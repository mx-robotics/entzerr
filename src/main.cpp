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
    GdkPixbuf *crop_pixbuf;
    GdkPixbuf *crop_scaled_pixbuf;
    double     crop_scale;
    int        crop_scaled_w, crop_scaled_h;
    // Shared
    GtkWidget *status_label;
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

static void update_crop(AppWidgets *app) {
    GdkPixbuf *warped = perspective_warp(app->original_pixbuf, app->pts);
    if (!warped) return;
    if (app->crop_pixbuf) g_object_unref(app->crop_pixbuf);
    app->crop_pixbuf = warped;
    update_crop_display(app);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(app->notebook), 1);
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
    app.crop_scroll = gtk_scrolled_window_new(nullptr, nullptr);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(app.crop_scroll),
                                   GTK_POLICY_NEVER, GTK_POLICY_NEVER);
    g_signal_connect(app.crop_scroll, "size-allocate",
                     G_CALLBACK(on_crop_scroll_size_allocate), &app);
    app.crop_drawing_area = gtk_drawing_area_new();
    g_signal_connect(app.crop_drawing_area, "draw", G_CALLBACK(on_crop_draw), &app);
    gtk_container_add(GTK_CONTAINER(app.crop_scroll), app.crop_drawing_area);
    gtk_notebook_append_page(GTK_NOTEBOOK(app.notebook), app.crop_scroll,
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

    if (app.original_pixbuf)   g_object_unref(app.original_pixbuf);
    if (app.scaled_pixbuf)     g_object_unref(app.scaled_pixbuf);
    if (app.crop_pixbuf)       g_object_unref(app.crop_pixbuf);
    if (app.crop_scaled_pixbuf) g_object_unref(app.crop_scaled_pixbuf);

    return 0;
}
