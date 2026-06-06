#include <gtk/gtk.h>
#include <cmath>
#include <cstring>

static constexpr int MAX_POINTS = 4;

struct AppWidgets {
    GtkWidget *window;
    GtkWidget *notebook;
    // Source tab
    GtkWidget *drawing_area;
    GtkWidget *scroll;
    GdkPixbuf *original_pixbuf;
    GdkPixbuf *scaled_pixbuf;
    double     scale;
    int        scaled_w, scaled_h;
    // Crop tab
    GtkWidget *crop_drawing_area;
    GtkWidget *crop_scroll;
    GtkWidget *combo_resize;
    GdkPixbuf *warp_pixbuf;         // raw perspective warp, never modified
    GdkPixbuf *clean_crop_pixbuf;   // warp + rotation + resize, no levels
    GdkPixbuf *crop_pixbuf;         // clean + levels applied (display + save)
    GdkPixbuf *crop_scaled_pixbuf;
    double     crop_scale;
    int        crop_scaled_w, crop_scaled_h;
    int        crop_rotation;       // 0 / 90 / 180 / 270
    int        resize_width;        // 0 = original
    bool       grayscale;
    // Histogram + levels
    GtkWidget *hist_area;
    int        hist_r[256];
    int        hist_g[256];
    int        hist_b[256];
    GtkWidget *scale_low;
    GtkWidget *scale_mid;
    GtkWidget *scale_high;
    int        level_low;           // input black point  0-254
    int        level_mid;           // input mid point    1-254
    int        level_high;          // input white point  1-255
    // Shared
    GtkWidget *status_label;
    char      *current_filename;
    double     pts[MAX_POINTS][2];   // original image coordinates
    int        point_count;
};

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

static void center_offset(int alloc_w, int alloc_h, int img_w, int img_h,
                           double &ox, double &oy) {
    ox = (alloc_w - img_w) / 2.0;
    oy = (alloc_h - img_h) / 2.0;
}

// Gaussian elimination with partial pivoting — solves 8×8 system Ax=b
static bool solve8(double A[8][8], double b[8], double x[8]) {
    double M[8][9];
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) M[i][j] = A[i][j];
        M[i][8] = b[i];
    }
    for (int col = 0; col < 8; col++) {
        int pivot = col;
        for (int row = col + 1; row < 8; row++)
            if (fabs(M[row][col]) > fabs(M[pivot][col])) pivot = row;
        if (pivot != col)
            for (int j = 0; j <= 8; j++) {
                double t = M[col][j]; M[col][j] = M[pivot][j]; M[pivot][j] = t;
            }
        if (fabs(M[col][col]) < 1e-10) return false;
        for (int row = 0; row < 8; row++) {
            if (row == col) continue;
            double f = M[row][col] / M[col][col];
            for (int j = col; j <= 8; j++) M[row][j] -= f * M[col][j];
        }
    }
    for (int i = 0; i < 8; i++) x[i] = M[i][8] / M[i][i];
    return true;
}

// Sort 4 arbitrary points into TL → TR → BR → BL order
static void sort_quad(double in[4][2], double out[4][2]) {
    double tmp[4][2];
    memcpy(tmp, in, sizeof(tmp));
    // Insertion sort by y
    for (int i = 1; i < 4; i++)
        for (int j = i; j > 0 && tmp[j][1] < tmp[j-1][1]; j--) {
            double tx = tmp[j][0], ty = tmp[j][1];
            tmp[j][0] = tmp[j-1][0]; tmp[j][1] = tmp[j-1][1];
            tmp[j-1][0] = tx; tmp[j-1][1] = ty;
        }
    // Top two: sort by x → TL, TR
    if (tmp[0][0] > tmp[1][0]) {
        double tx = tmp[0][0], ty = tmp[0][1];
        tmp[0][0] = tmp[1][0]; tmp[0][1] = tmp[1][1];
        tmp[1][0] = tx; tmp[1][1] = ty;
    }
    // Bottom two: sort by x → BL, BR
    if (tmp[2][0] > tmp[3][0]) {
        double tx = tmp[2][0], ty = tmp[2][1];
        tmp[2][0] = tmp[3][0]; tmp[2][1] = tmp[3][1];
        tmp[3][0] = tx; tmp[3][1] = ty;
    }
    memcpy(out[0], tmp[0], 2 * sizeof(double)); // TL
    memcpy(out[1], tmp[1], 2 * sizeof(double)); // TR
    memcpy(out[2], tmp[3], 2 * sizeof(double)); // BR
    memcpy(out[3], tmp[2], 2 * sizeof(double)); // BL
}

// Perspective warp: extract the quadrilateral defined by pts into a rectangle
static GdkPixbuf *perspective_warp(GdkPixbuf *src, double pts[4][2]) {
    double sp[4][2];
    sort_quad(pts, sp); // TL, TR, BR, BL

    // Output size = average of opposite edge lengths
    double w = (hypot(sp[1][0]-sp[0][0], sp[1][1]-sp[0][1]) +
                hypot(sp[2][0]-sp[3][0], sp[2][1]-sp[3][1])) / 2.0;
    double h = (hypot(sp[3][0]-sp[0][0], sp[3][1]-sp[0][1]) +
                hypot(sp[2][0]-sp[1][0], sp[2][1]-sp[1][1])) / 2.0;
    int dst_w = MAX(1, (int)w);
    int dst_h = MAX(1, (int)h);

    // Destination corners: TL, TR, BR, BL
    double dp[4][2] = {
        {0, 0}, {(double)(dst_w-1), 0},
        {(double)(dst_w-1), (double)(dst_h-1)}, {0, (double)(dst_h-1)}
    };

    // Build inverse homography H: dst → src
    // sx = (h0*dx + h1*dy + h2) / (h6*dx + h7*dy + 1)
    // sy = (h3*dx + h4*dy + h5) / (h6*dx + h7*dy + 1)
    double A[8][8] = {}, b[8] = {};
    for (int i = 0; i < 4; i++) {
        double dx = dp[i][0], dy = dp[i][1];
        double sx = sp[i][0], sy = sp[i][1];
        A[2*i][0]=dx; A[2*i][1]=dy; A[2*i][2]=1;
        A[2*i][6]=-sx*dx; A[2*i][7]=-sx*dy; b[2*i]=sx;
        A[2*i+1][3]=dx; A[2*i+1][4]=dy; A[2*i+1][5]=1;
        A[2*i+1][6]=-sy*dx; A[2*i+1][7]=-sy*dy; b[2*i+1]=sy;
    }
    double hh[8];
    if (!solve8(A, b, hh)) return nullptr;

    GdkPixbuf *dst = gdk_pixbuf_new(GDK_COLORSPACE_RGB,
                                    gdk_pixbuf_get_has_alpha(src), 8, dst_w, dst_h);
    if (!dst) return nullptr;

    int src_w      = gdk_pixbuf_get_width(src);
    int src_h      = gdk_pixbuf_get_height(src);
    int src_ch     = gdk_pixbuf_get_n_channels(src);
    int src_stride = gdk_pixbuf_get_rowstride(src);
    int dst_stride = gdk_pixbuf_get_rowstride(dst);
    const guchar *sp_px = gdk_pixbuf_get_pixels(src);
    guchar       *dp_px = gdk_pixbuf_get_pixels(dst);

    for (int dy_i = 0; dy_i < dst_h; dy_i++) {
        guchar *drow = dp_px + dy_i * dst_stride;
        for (int dx_i = 0; dx_i < dst_w; dx_i++) {
            double ww = hh[6]*dx_i + hh[7]*dy_i + 1.0;
            double sx  = (hh[0]*dx_i + hh[1]*dy_i + hh[2]) / ww;
            double sy  = (hh[3]*dx_i + hh[4]*dy_i + hh[5]) / ww;
            guchar *out = drow + dx_i * src_ch;
            int x0 = (int)sx, y0 = (int)sy;
            if (x0 < 0 || y0 < 0 || x0+1 >= src_w || y0+1 >= src_h) {
                for (int c = 0; c < src_ch; c++) out[c] = 0;
                continue;
            }
            double fx = sx - x0, fy = sy - y0;
            const guchar *r00 = sp_px + y0*src_stride     + x0*src_ch;
            const guchar *r10 = sp_px + y0*src_stride     + (x0+1)*src_ch;
            const guchar *r01 = sp_px + (y0+1)*src_stride + x0*src_ch;
            const guchar *r11 = sp_px + (y0+1)*src_stride + (x0+1)*src_ch;
            for (int c = 0; c < src_ch; c++)
                out[c] = (guchar)(r00[c]*(1-fx)*(1-fy) + r10[c]*fx*(1-fy) +
                                  r01[c]*(1-fx)*fy     + r11[c]*fx*fy + 0.5);
        }
    }
    return dst;
}

// ---------------------------------------------------------------------------
// Source tab — scale & draw
// ---------------------------------------------------------------------------

static void update_scale(AppWidgets *app) {
    if (!app->original_pixbuf) return;
    int area_w = gtk_widget_get_allocated_width(app->scroll);
    int area_h = gtk_widget_get_allocated_height(app->scroll);
    if (area_w <= 1 || area_h <= 1) return;
    int img_w = gdk_pixbuf_get_width(app->original_pixbuf);
    int img_h = gdk_pixbuf_get_height(app->original_pixbuf);
    app->scale    = MIN((double)area_w / img_w, (double)area_h / img_h);
    app->scaled_w = MAX(1, (int)(img_w * app->scale));
    app->scaled_h = MAX(1, (int)(img_h * app->scale));
    if (app->scaled_pixbuf) g_object_unref(app->scaled_pixbuf);
    app->scaled_pixbuf = gdk_pixbuf_scale_simple(
        app->original_pixbuf, app->scaled_w, app->scaled_h, GDK_INTERP_BILINEAR);
    gtk_widget_queue_draw(app->drawing_area);
}

static void on_scroll_size_allocate(GtkWidget *, GdkRectangle *, gpointer data) {
    update_scale(static_cast<AppWidgets *>(data));
}

static gboolean on_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    if (!app->scaled_pixbuf) return FALSE;

    double ox, oy;
    center_offset(gtk_widget_get_allocated_width(widget),
                  gtk_widget_get_allocated_height(widget),
                  app->scaled_w, app->scaled_h, ox, oy);

    gdk_cairo_set_source_pixbuf(cr, app->scaled_pixbuf, ox, oy);
    cairo_paint(cr);

    if (app->point_count == 0) return FALSE;

    double dx[MAX_POINTS], dy[MAX_POINTS];
    for (int i = 0; i < app->point_count; i++) {
        dx[i] = ox + app->pts[i][0] * app->scale;
        dy[i] = oy + app->pts[i][1] * app->scale;
    }

    cairo_set_source_rgba(cr, 0.0, 0.85, 1.0, 0.9);
    cairo_set_line_width(cr, 2.0);
    cairo_move_to(cr, dx[0], dy[0]);
    for (int i = 1; i < app->point_count; i++)
        cairo_line_to(cr, dx[i], dy[i]);
    if (app->point_count == MAX_POINTS)
        cairo_close_path(cr);
    cairo_stroke(cr);

    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 12.0);
    for (int i = 0; i < app->point_count; i++) {
        cairo_set_source_rgba(cr, 1.0, 0.3, 0.0, 1.0);
        cairo_arc(cr, dx[i], dy[i], 5.0, 0, 2 * G_PI);
        cairo_fill(cr);
        char label[4];
        g_snprintf(label, sizeof(label), "%d", i + 1);
        cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);
        cairo_move_to(cr, dx[i] + 8, dy[i] + 5);
        cairo_show_text(cr, label);
    }
    return FALSE;
}

// ---------------------------------------------------------------------------
// Crop tab — warp, scale & draw
// ---------------------------------------------------------------------------

static void update_crop_display(AppWidgets *app) {
    if (!app->crop_pixbuf) return;
    int area_w = gtk_widget_get_allocated_width(app->crop_scroll);
    int area_h = gtk_widget_get_allocated_height(app->crop_scroll);
    if (area_w <= 1 || area_h <= 1) return;
    int img_w = gdk_pixbuf_get_width(app->crop_pixbuf);
    int img_h = gdk_pixbuf_get_height(app->crop_pixbuf);
    app->crop_scale    = MIN((double)area_w / img_w, (double)area_h / img_h);
    app->crop_scaled_w = MAX(1, (int)(img_w * app->crop_scale));
    app->crop_scaled_h = MAX(1, (int)(img_h * app->crop_scale));
    if (app->crop_scaled_pixbuf) g_object_unref(app->crop_scaled_pixbuf);
    app->crop_scaled_pixbuf = gdk_pixbuf_scale_simple(
        app->crop_pixbuf, app->crop_scaled_w, app->crop_scaled_h, GDK_INTERP_BILINEAR);
    gtk_widget_queue_draw(app->crop_drawing_area);
}

static void on_crop_scroll_size_allocate(GtkWidget *, GdkRectangle *, gpointer data) {
    update_crop_display(static_cast<AppWidgets *>(data));
}

static gboolean on_crop_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    if (!app->crop_scaled_pixbuf) return FALSE;
    double ox, oy;
    center_offset(gtk_widget_get_allocated_width(widget),
                  gtk_widget_get_allocated_height(widget),
                  app->crop_scaled_w, app->crop_scaled_h, ox, oy);
    gdk_cairo_set_source_pixbuf(cr, app->crop_scaled_pixbuf, ox, oy);
    cairo_paint(cr);
    return FALSE;
}

static GdkPixbuf *apply_levels(GdkPixbuf *src, int low, int mid, int high) {
    low  = MAX(0,       MIN(253, low));
    high = MAX(low + 2, MIN(255, high));
    mid  = MAX(low + 1, MIN(high - 1, mid));

    double mid_frac = (double)(mid - low) / (high - low);
    mid_frac = MAX(1e-4, MIN(1.0 - 1e-4, mid_frac));
    double gamma = log(0.5) / log(mid_frac);

    guchar lut[256];
    for (int i = 0; i < 256; i++) {
        double v = (i - low) / (double)(high - low);
        v = MAX(0.0, MIN(1.0, v));
        lut[i] = (guchar)(pow(v, gamma) * 255.0 + 0.5);
    }

    int w      = gdk_pixbuf_get_width(src);
    int h      = gdk_pixbuf_get_height(src);
    int ch     = gdk_pixbuf_get_n_channels(src);
    int stride = gdk_pixbuf_get_rowstride(src);
    GdkPixbuf *dst = gdk_pixbuf_copy(src);
    guchar *px = gdk_pixbuf_get_pixels(dst);
    for (int y = 0; y < h; y++) {
        guchar *row = px + y * stride;
        for (int x = 0; x < w; x++)
            for (int c = 0; c < MIN(ch, 3); c++)
                row[x * ch + c] = lut[row[x * ch + c]];
    }
    return dst;
}

static GdkPixbuf *apply_grayscale(GdkPixbuf *src) {
    GdkPixbuf *dst = gdk_pixbuf_copy(src);
    int w      = gdk_pixbuf_get_width(dst);
    int h      = gdk_pixbuf_get_height(dst);
    int ch     = gdk_pixbuf_get_n_channels(dst);
    int stride = gdk_pixbuf_get_rowstride(dst);
    guchar *px = gdk_pixbuf_get_pixels(dst);
    for (int y = 0; y < h; y++) {
        guchar *row = px + y * stride;
        for (int x = 0; x < w; x++) {
            guchar *p = row + x * ch;
            guchar lum = (guchar)(0.299 * p[0] + 0.587 * p[1] + 0.114 * p[2] + 0.5);
            p[0] = p[1] = p[2] = lum;
        }
    }
    return dst;
}

static void compute_histogram(AppWidgets *app) {
    memset(app->hist_r, 0, sizeof(app->hist_r));
    memset(app->hist_g, 0, sizeof(app->hist_g));
    memset(app->hist_b, 0, sizeof(app->hist_b));
    if (!app->clean_crop_pixbuf) return;

    int ch     = gdk_pixbuf_get_n_channels(app->clean_crop_pixbuf);
    int stride = gdk_pixbuf_get_rowstride(app->clean_crop_pixbuf);
    int w      = gdk_pixbuf_get_width(app->clean_crop_pixbuf);
    int h      = gdk_pixbuf_get_height(app->clean_crop_pixbuf);
    const guchar *px = gdk_pixbuf_get_pixels(app->clean_crop_pixbuf);

    for (int y = 0; y < h; y++) {
        const guchar *row = px + y * stride;
        for (int x = 0; x < w; x++) {
            app->hist_r[row[x * ch + 0]]++;
            app->hist_g[row[x * ch + 1]]++;
            app->hist_b[row[x * ch + 2]]++;
        }
    }
}

static gboolean on_hist_draw(GtkWidget *widget, cairo_t *cr, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);

    int aw = gtk_widget_get_allocated_width(widget);
    int ah = gtk_widget_get_allocated_height(widget);

    // Background
    cairo_set_source_rgb(cr, 0.12, 0.12, 0.12);
    cairo_paint(cr);

    if (!app->crop_pixbuf) return FALSE;

    // Find global max for normalisation
    int max_val = 1;
    for (int i = 0; i < 256; i++) {
        max_val = MAX(max_val, app->hist_r[i]);
        max_val = MAX(max_val, app->hist_g[i]);
        max_val = MAX(max_val, app->hist_b[i]);
    }

    const int pad = 4;
    double plot_w = aw - 2 * pad;
    double plot_h = ah - 2 * pad;
    double bin_w  = plot_w / 256.0;

    struct { int *data; double r, g, b; } channels[3] = {
        { app->hist_r, 0.9, 0.2, 0.2 },
        { app->hist_g, 0.2, 0.85, 0.2 },
        { app->hist_b, 0.3, 0.55, 1.0 },
    };

    for (auto &ch : channels) {
        cairo_set_source_rgba(cr, ch.r, ch.g, ch.b, 0.55);
        for (int i = 0; i < 256; i++) {
            double bar_h = (double)ch.data[i] / max_val * plot_h;
            double x     = pad + i * bin_w;
            cairo_rectangle(cr, x, pad + plot_h - bar_h, MAX(1.0, bin_w), bar_h);
        }
        cairo_fill(cr);
    }

    // Level marker lines
    struct { int val; double r, g, b; } markers[3] = {
        { app->level_low,  1.0, 1.0, 0.0 },
        { app->level_mid,  0.8, 0.8, 0.8 },
        { app->level_high, 1.0, 1.0, 0.0 },
    };
    cairo_set_line_width(cr, 1.5);
    for (auto &m : markers) {
        double x = pad + (double)m.val / 255.0 * plot_w;
        cairo_set_source_rgba(cr, m.r, m.g, m.b, 0.9);
        cairo_move_to(cr, x, pad);
        cairo_line_to(cr, x, pad + plot_h);
        cairo_stroke(cr);
    }

    // Border
    cairo_set_source_rgba(cr, 0.4, 0.4, 0.4, 1.0);
    cairo_set_line_width(cr, 1.0);
    cairo_rectangle(cr, pad, pad, plot_w, plot_h);
    cairo_stroke(cr);

    return FALSE;
}

static void recompute_crop(AppWidgets *app) {
    if (!app->warp_pixbuf) return;

    GdkPixbuf *tmp = app->warp_pixbuf;
    g_object_ref(tmp);

    if (app->crop_rotation != 0) {
        GdkPixbufRotation rot;
        switch (app->crop_rotation) {
            case  90: rot = GDK_PIXBUF_ROTATE_CLOCKWISE;        break;
            case 180: rot = GDK_PIXBUF_ROTATE_UPSIDEDOWN;       break;
            default:  rot = GDK_PIXBUF_ROTATE_COUNTERCLOCKWISE; break;
        }
        GdkPixbuf *rotated = gdk_pixbuf_rotate_simple(tmp, rot);
        g_object_unref(tmp);
        tmp = rotated;
    }

    if (app->resize_width > 0 && tmp) {
        int w = gdk_pixbuf_get_width(tmp);
        int h = gdk_pixbuf_get_height(tmp);
        int new_h = MAX(1, (int)((double)h * app->resize_width / w));
        GdkPixbuf *resized = gdk_pixbuf_scale_simple(tmp, app->resize_width, new_h,
                                                      GDK_INTERP_BILINEAR);
        g_object_unref(tmp);
        tmp = resized;
    }

    if (app->grayscale) {
        GdkPixbuf *gray = apply_grayscale(tmp);
        g_object_unref(tmp);
        tmp = gray;
    }

    if (app->clean_crop_pixbuf) g_object_unref(app->clean_crop_pixbuf);
    app->clean_crop_pixbuf = tmp;

    GdkPixbuf *leveled = apply_levels(tmp, app->level_low, app->level_mid, app->level_high);
    if (app->crop_pixbuf) g_object_unref(app->crop_pixbuf);
    app->crop_pixbuf = leveled;

    update_crop_display(app);
    compute_histogram(app);
    if (app->hist_area) gtk_widget_queue_draw(app->hist_area);
}

static void on_grayscale_toggled(GtkToggleButton *btn, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    app->grayscale = gtk_toggle_button_get_active(btn);
    recompute_crop(app);
}

static void on_level_changed(GtkRange *, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    app->level_low  = (int)gtk_range_get_value(GTK_RANGE(app->scale_low));
    app->level_mid  = (int)gtk_range_get_value(GTK_RANGE(app->scale_mid));
    app->level_high = (int)gtk_range_get_value(GTK_RANGE(app->scale_high));
    recompute_crop(app);
}

static void update_crop(AppWidgets *app) {
    GdkPixbuf *warped = perspective_warp(app->original_pixbuf, app->pts);
    if (!warped) return;
    if (app->warp_pixbuf) g_object_unref(app->warp_pixbuf);
    app->warp_pixbuf   = warped;
    app->crop_rotation = 0;
    app->level_low     = 0;
    app->level_mid     = 128;
    app->level_high    = 255;
    if (app->scale_low) {
        g_signal_handlers_block_by_func(app->scale_low,  (gpointer)on_level_changed, app);
        g_signal_handlers_block_by_func(app->scale_mid,  (gpointer)on_level_changed, app);
        g_signal_handlers_block_by_func(app->scale_high, (gpointer)on_level_changed, app);
        gtk_range_set_value(GTK_RANGE(app->scale_low),  0);
        gtk_range_set_value(GTK_RANGE(app->scale_mid),  128);
        gtk_range_set_value(GTK_RANGE(app->scale_high), 255);
        g_signal_handlers_unblock_by_func(app->scale_low,  (gpointer)on_level_changed, app);
        g_signal_handlers_unblock_by_func(app->scale_mid,  (gpointer)on_level_changed, app);
        g_signal_handlers_unblock_by_func(app->scale_high, (gpointer)on_level_changed, app);
    }
    recompute_crop(app);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(app->notebook), 1);
}

static void on_rotate_left(GtkButton *, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    app->crop_rotation = (app->crop_rotation + 270) % 360;
    recompute_crop(app);
}

static void on_rotate_right(GtkButton *, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    app->crop_rotation = (app->crop_rotation + 90) % 360;
    recompute_crop(app);
}

static void on_resize_changed(GtkComboBoxText *combo, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    char *text = gtk_combo_box_text_get_active_text(combo);
    if (!text) return;
    app->resize_width = (g_strcmp0(text, "Original") == 0) ? 0 : atoi(text);
    g_free(text);
    recompute_crop(app);
}

// ---------------------------------------------------------------------------
// Save helpers
// ---------------------------------------------------------------------------

static const char *format_from_filename(const char *filename) {
    const char *ext = strrchr(filename, '.');
    if (!ext) return "png";
    ext++;
    if (g_ascii_strcasecmp(ext, "jpg") == 0 || g_ascii_strcasecmp(ext, "jpeg") == 0) return "jpeg";
    if (g_ascii_strcasecmp(ext, "png")  == 0) return "png";
    if (g_ascii_strcasecmp(ext, "bmp")  == 0) return "bmp";
    if (g_ascii_strcasecmp(ext, "tiff") == 0 || g_ascii_strcasecmp(ext, "tif") == 0) return "tiff";
    return "png";
}

static char *make_suffixed_filename(const char *filename, const char *suffix) {
    const char *dot = strrchr(filename, '.');
    if (!dot) return g_strdup_printf("%s%s", filename, suffix);
    char *base = g_strndup(filename, dot - filename);
    char *result = g_strdup_printf("%s%s%s", base, suffix, dot);
    g_free(base);
    return result;
}

static char *make_crop_filename(const char *filename) {
    return make_suffixed_filename(filename, "_crop");
}

static void save_pixbuf_to(AppWidgets *app, GdkPixbuf *pixbuf, const char *filename) {
    const char *fmt = format_from_filename(filename);
    GError *error = nullptr;
    gboolean ok;
    if (g_strcmp0(fmt, "jpeg") == 0)
        ok = gdk_pixbuf_save(pixbuf, filename, fmt, &error, "quality", "95", nullptr);
    else
        ok = gdk_pixbuf_save(pixbuf, filename, fmt, &error, nullptr);

    if (ok) {
        char buf[256];
        g_snprintf(buf, sizeof(buf), "Saved: %s", g_path_get_basename(filename));
        gtk_label_set_text(GTK_LABEL(app->status_label), buf);
    } else {
        GtkWidget *dlg = gtk_message_dialog_new(GTK_WINDOW(app->window),
            GTK_DIALOG_MODAL, GTK_MESSAGE_ERROR, GTK_BUTTONS_OK,
            "Save failed: %s", error ? error->message : "unknown error");
        gtk_dialog_run(GTK_DIALOG(dlg));
        gtk_widget_destroy(dlg);
        if (error) g_error_free(error);
    }
}

static void on_save_replace(GtkButton *, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    if (!app->crop_pixbuf || !app->current_filename) return;
    char *org_filename = make_suffixed_filename(app->current_filename, "-org");
    save_pixbuf_to(app, app->original_pixbuf, org_filename);
    g_free(org_filename);
    save_pixbuf_to(app, app->crop_pixbuf, app->current_filename);
}

static void on_save_as(GtkButton *, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    if (!app->crop_pixbuf) return;

    GtkWidget *dialog = gtk_file_chooser_dialog_new(
        "Save Cropped Image", GTK_WINDOW(app->window), GTK_FILE_CHOOSER_ACTION_SAVE,
        "_Cancel", GTK_RESPONSE_CANCEL, "_Save", GTK_RESPONSE_ACCEPT, nullptr);
    gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);

    if (app->current_filename) {
        char *suggested = make_crop_filename(app->current_filename);
        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog),
                                            g_path_get_dirname(app->current_filename));
        gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog),
                                          g_path_get_basename(suggested));
        g_free(suggested);
    }

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        save_pixbuf_to(app, app->crop_pixbuf, filename);
        g_free(filename);
    }
    gtk_widget_destroy(dialog);
}

// ---------------------------------------------------------------------------
// Click handler
// ---------------------------------------------------------------------------

static gboolean on_click(GtkWidget *widget, GdkEventButton *event, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    if (!app->original_pixbuf || event->button != 1) return FALSE;

    double ox, oy;
    center_offset(gtk_widget_get_allocated_width(widget),
                  gtk_widget_get_allocated_height(widget),
                  app->scaled_w, app->scaled_h, ox, oy);

    double img_x = (event->x - ox) / app->scale;
    double img_y = (event->y - oy) / app->scale;

    int orig_w = gdk_pixbuf_get_width(app->original_pixbuf);
    int orig_h = gdk_pixbuf_get_height(app->original_pixbuf);
    if (img_x < 0 || img_y < 0 || img_x >= orig_w || img_y >= orig_h) {
        gtk_label_set_text(GTK_LABEL(app->status_label), "Outside image");
        return TRUE;
    }

    if (app->point_count == MAX_POINTS) app->point_count = 0;

    app->pts[app->point_count][0] = img_x;
    app->pts[app->point_count][1] = img_y;
    app->point_count++;

    char buf[128];
    if (app->point_count < MAX_POINTS)
        g_snprintf(buf, sizeof(buf), "Point %d: x=%.0f  y=%.0f  (%d / %d)",
                   app->point_count, img_x, img_y, app->point_count, MAX_POINTS);
    else
        g_snprintf(buf, sizeof(buf), "Point %d: x=%.0f  y=%.0f  — generating crop…",
                   app->point_count, img_x, img_y);

    gtk_label_set_text(GTK_LABEL(app->status_label), buf);
    gtk_widget_queue_draw(app->drawing_area);

    if (app->point_count == MAX_POINTS)
        update_crop(app);

    return TRUE;
}

// ---------------------------------------------------------------------------
// Image loading
// ---------------------------------------------------------------------------

static void load_image(AppWidgets *app, const char *filename) {
    GError *error = nullptr;
    GdkPixbuf *pixbuf = gdk_pixbuf_new_from_file(filename, &error);
    if (pixbuf) {
        if (app->original_pixbuf) g_object_unref(app->original_pixbuf);
        app->original_pixbuf = pixbuf;
        app->point_count = 0;
        g_free(app->current_filename);
        app->current_filename = g_strdup(filename);
        gtk_notebook_set_current_page(GTK_NOTEBOOK(app->notebook), 0);
        update_scale(app);
        gtk_window_set_title(GTK_WINDOW(app->window), g_path_get_basename(filename));
        gtk_label_set_text(GTK_LABEL(app->status_label),
                           "Click 4 points on the image to extract a perspective-corrected crop");
    } else {
        g_warning("Could not load image: %s", error->message);
        g_error_free(error);
    }
}

static void on_open(GtkMenuItem *, gpointer data) {
    AppWidgets *app = static_cast<AppWidgets *>(data);
    GtkWidget *dialog = gtk_file_chooser_dialog_new(
        "Open Image", GTK_WINDOW(app->window), GTK_FILE_CHOOSER_ACTION_OPEN,
        "_Cancel", GTK_RESPONSE_CANCEL, "_Open", GTK_RESPONSE_ACCEPT, nullptr);

    GtkFileFilter *img_filter = gtk_file_filter_new();
    gtk_file_filter_set_name(img_filter, "Images");
    gtk_file_filter_add_mime_type(img_filter, "image/png");
    gtk_file_filter_add_mime_type(img_filter, "image/jpeg");
    gtk_file_filter_add_mime_type(img_filter, "image/gif");
    gtk_file_filter_add_mime_type(img_filter, "image/bmp");
    gtk_file_filter_add_mime_type(img_filter, "image/webp");
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), img_filter);

    GtkFileFilter *all_filter = gtk_file_filter_new();
    gtk_file_filter_set_name(all_filter, "All files");
    gtk_file_filter_add_pattern(all_filter, "*");
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), all_filter);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT) {
        char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        load_image(app, filename);
        g_free(filename);
    }
    gtk_widget_destroy(dialog);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    gtk_init(&argc, &argv);

    AppWidgets app = {};

    app.window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(app.window), "Image Viewer");
    gtk_window_set_default_size(GTK_WINDOW(app.window), 800, 600);
    g_signal_connect(app.window, "destroy", G_CALLBACK(gtk_main_quit), nullptr);

    GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_container_add(GTK_CONTAINER(app.window), vbox);

    // Menu bar
    GtkWidget *menubar   = gtk_menu_bar_new();
    GtkWidget *file_menu = gtk_menu_new();
    GtkWidget *file_item = gtk_menu_item_new_with_mnemonic("_File");
    gtk_menu_item_set_submenu(GTK_MENU_ITEM(file_item), file_menu);
    gtk_menu_shell_append(GTK_MENU_SHELL(menubar), file_item);

    GtkWidget *open_item = gtk_menu_item_new_with_mnemonic("_Open...");
    g_signal_connect(open_item, "activate", G_CALLBACK(on_open), &app);
    gtk_menu_shell_append(GTK_MENU_SHELL(file_menu), open_item);
    gtk_menu_shell_append(GTK_MENU_SHELL(file_menu), gtk_separator_menu_item_new());
    GtkWidget *quit_item = gtk_menu_item_new_with_mnemonic("_Quit");
    g_signal_connect(quit_item, "activate", G_CALLBACK(gtk_main_quit), nullptr);
    gtk_menu_shell_append(GTK_MENU_SHELL(file_menu), quit_item);
    gtk_box_pack_start(GTK_BOX(vbox), menubar, FALSE, FALSE, 0);

    // Notebook
    app.notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(vbox), app.notebook, TRUE, TRUE, 0);

    // Tab 1 — Source
    app.scroll = gtk_scrolled_window_new(nullptr, nullptr);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(app.scroll),
                                   GTK_POLICY_NEVER, GTK_POLICY_NEVER);
    g_signal_connect(app.scroll, "size-allocate", G_CALLBACK(on_scroll_size_allocate), &app);
    app.drawing_area = gtk_drawing_area_new();
    gtk_widget_add_events(app.drawing_area, GDK_BUTTON_PRESS_MASK);
    g_signal_connect(app.drawing_area, "draw",               G_CALLBACK(on_draw),  &app);
    g_signal_connect(app.drawing_area, "button-press-event", G_CALLBACK(on_click), &app);
    gtk_container_add(GTK_CONTAINER(app.scroll), app.drawing_area);
    gtk_notebook_append_page(GTK_NOTEBOOK(app.notebook), app.scroll,
                             gtk_label_new("Source"));

    // Tab 2 — Crop
    GtkWidget *crop_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);

    GtkWidget *toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
    gtk_widget_set_margin_start(toolbar, 4);
    gtk_widget_set_margin_end(toolbar, 4);
    gtk_widget_set_margin_top(toolbar, 4);
    gtk_widget_set_margin_bottom(toolbar, 4);

    GtkWidget *btn_rotate_left  = gtk_button_new_with_label("↺  Rotate Left");
    GtkWidget *btn_rotate_right = gtk_button_new_with_label("↻  Rotate Right");
    g_signal_connect(btn_rotate_left,  "clicked", G_CALLBACK(on_rotate_left),  &app);
    g_signal_connect(btn_rotate_right, "clicked", G_CALLBACK(on_rotate_right), &app);
    gtk_box_pack_start(GTK_BOX(toolbar), btn_rotate_left,  FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(toolbar), btn_rotate_right, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(toolbar),
                       gtk_separator_new(GTK_ORIENTATION_VERTICAL), FALSE, FALSE, 4);

    GtkWidget *btn_gray = gtk_toggle_button_new_with_label("Grayscale");
    g_signal_connect(btn_gray, "toggled", G_CALLBACK(on_grayscale_toggled), &app);
    gtk_box_pack_start(GTK_BOX(toolbar), btn_gray, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(toolbar),
                       gtk_separator_new(GTK_ORIENTATION_VERTICAL), FALSE, FALSE, 4);

    gtk_box_pack_start(GTK_BOX(toolbar), gtk_label_new("Resize:"), FALSE, FALSE, 0);

    app.combo_resize = gtk_combo_box_text_new();
    const char *sizes[] = { "400 px", "500 px", "800 px", "1000 px", "1500 px", "Original", nullptr };
    for (int i = 0; sizes[i]; i++)
        gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(app.combo_resize), sizes[i]);
    gtk_combo_box_set_active(GTK_COMBO_BOX(app.combo_resize), 5); // Original
    g_signal_connect(app.combo_resize, "changed", G_CALLBACK(on_resize_changed), &app);
    gtk_box_pack_start(GTK_BOX(toolbar), app.combo_resize, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(crop_vbox), toolbar, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(crop_vbox),
                       gtk_separator_new(GTK_ORIENTATION_HORIZONTAL), FALSE, FALSE, 0);

    GtkWidget *content_pane = gtk_paned_new(GTK_ORIENTATION_HORIZONTAL);
    gtk_box_pack_start(GTK_BOX(crop_vbox), content_pane, TRUE, TRUE, 0);

    app.crop_scroll = gtk_scrolled_window_new(nullptr, nullptr);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(app.crop_scroll),
                                   GTK_POLICY_NEVER, GTK_POLICY_NEVER);
    g_signal_connect(app.crop_scroll, "size-allocate",
                     G_CALLBACK(on_crop_scroll_size_allocate), &app);
    app.crop_drawing_area = gtk_drawing_area_new();
    g_signal_connect(app.crop_drawing_area, "draw", G_CALLBACK(on_crop_draw), &app);
    gtk_container_add(GTK_CONTAINER(app.crop_scroll), app.crop_drawing_area);
    gtk_paned_pack1(GTK_PANED(content_pane), app.crop_scroll, TRUE, FALSE);

    GtkWidget *hist_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
    gtk_widget_set_margin_start(hist_vbox, 4);
    gtk_widget_set_margin_end(hist_vbox, 4);
    gtk_widget_set_margin_top(hist_vbox, 4);
    gtk_widget_set_margin_bottom(hist_vbox, 4);
    gtk_widget_set_size_request(hist_vbox, 220, -1);

    GtkWidget *hist_frame = gtk_frame_new("Histogram");
    app.hist_area = gtk_drawing_area_new();
    g_signal_connect(app.hist_area, "draw", G_CALLBACK(on_hist_draw), &app);
    gtk_container_add(GTK_CONTAINER(hist_frame), app.hist_area);
    gtk_box_pack_start(GTK_BOX(hist_vbox), hist_frame, TRUE, TRUE, 0);

    GtkWidget *levels_frame = gtk_frame_new("Levels");
    GtkWidget *levels_grid = gtk_grid_new();
    gtk_grid_set_column_spacing(GTK_GRID(levels_grid), 4);
    gtk_grid_set_row_spacing(GTK_GRID(levels_grid), 2);
    gtk_widget_set_margin_start(levels_grid, 4);
    gtk_widget_set_margin_end(levels_grid, 4);
    gtk_widget_set_margin_top(levels_grid, 4);
    gtk_widget_set_margin_bottom(levels_grid, 4);

    struct { const char *label; GtkWidget **scale; double val; } level_rows[3] = {
        { "Low",  &app.scale_low,  0   },
        { "Mid",  &app.scale_mid,  128 },
        { "High", &app.scale_high, 255 },
    };
    for (int i = 0; i < 3; i++) {
        GtkWidget *lbl = gtk_label_new(level_rows[i].label);
        gtk_label_set_xalign(GTK_LABEL(lbl), 1.0);
        gtk_grid_attach(GTK_GRID(levels_grid), lbl, 0, i, 1, 1);

        *level_rows[i].scale = gtk_scale_new_with_range(
            GTK_ORIENTATION_HORIZONTAL, 0, 255, 1);
        gtk_scale_set_draw_value(GTK_SCALE(*level_rows[i].scale), TRUE);
        gtk_scale_set_value_pos(GTK_SCALE(*level_rows[i].scale), GTK_POS_RIGHT);
        gtk_range_set_value(GTK_RANGE(*level_rows[i].scale), level_rows[i].val);
        gtk_widget_set_hexpand(*level_rows[i].scale, TRUE);
        g_signal_connect(*level_rows[i].scale, "value-changed",
                         G_CALLBACK(on_level_changed), &app);
        gtk_grid_attach(GTK_GRID(levels_grid), *level_rows[i].scale, 1, i, 1, 1);
    }
    gtk_container_add(GTK_CONTAINER(levels_frame), levels_grid);
    gtk_box_pack_start(GTK_BOX(hist_vbox), levels_frame, FALSE, FALSE, 0);

    gtk_paned_pack2(GTK_PANED(content_pane), hist_vbox, FALSE, FALSE);

    GtkWidget *savebar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
    gtk_widget_set_margin_start(savebar, 4);
    gtk_widget_set_margin_end(savebar, 4);
    gtk_widget_set_margin_top(savebar, 4);
    gtk_widget_set_margin_bottom(savebar, 4);
    GtkWidget *btn_save_replace = gtk_button_new_with_label("Replace");
    GtkWidget *btn_save_as      = gtk_button_new_with_label("Save as");
    g_signal_connect(btn_save_replace, "clicked", G_CALLBACK(on_save_replace), &app);
    g_signal_connect(btn_save_as,      "clicked", G_CALLBACK(on_save_as),      &app);
    gtk_box_pack_end(GTK_BOX(savebar), btn_save_as,      FALSE, FALSE, 0);
    gtk_box_pack_end(GTK_BOX(savebar), btn_save_replace, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(crop_vbox),
                       gtk_separator_new(GTK_ORIENTATION_HORIZONTAL), FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(crop_vbox), savebar, FALSE, FALSE, 0);

    gtk_notebook_append_page(GTK_NOTEBOOK(app.notebook), crop_vbox,
                             gtk_label_new("Crop"));

    // Status bar
    gtk_box_pack_start(GTK_BOX(vbox),
                       gtk_separator_new(GTK_ORIENTATION_HORIZONTAL), FALSE, FALSE, 0);
    GtkWidget *statusbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
    gtk_box_pack_start(GTK_BOX(vbox), statusbar, FALSE, FALSE, 0);
    app.status_label = gtk_label_new("Open an image to get started");
    gtk_label_set_xalign(GTK_LABEL(app.status_label), 0.0);
    gtk_widget_set_margin_start(app.status_label, 6);
    gtk_widget_set_margin_top(app.status_label, 3);
    gtk_widget_set_margin_bottom(app.status_label, 3);
    gtk_box_pack_start(GTK_BOX(statusbar), app.status_label, TRUE, TRUE, 0);

    gtk_widget_show_all(app.window);

    if (argc == 2)
        load_image(&app, argv[1]);

    gtk_main();

    if (app.original_pixbuf)    g_object_unref(app.original_pixbuf);
    if (app.scaled_pixbuf)      g_object_unref(app.scaled_pixbuf);
    if (app.warp_pixbuf)        g_object_unref(app.warp_pixbuf);
    if (app.clean_crop_pixbuf)  g_object_unref(app.clean_crop_pixbuf);
    if (app.crop_pixbuf)        g_object_unref(app.crop_pixbuf);
    if (app.crop_scaled_pixbuf) g_object_unref(app.crop_scaled_pixbuf);
    g_free(app.current_filename);

    return 0;
}
