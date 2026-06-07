#include <entzerr/entzerr.hpp>
#include <cmath>
#include <cstring>

void center_offset(int alloc_w, int alloc_h, int img_w, int img_h,
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
GdkPixbuf *perspective_warp(GdkPixbuf *src, double pts[4][2]) {
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
