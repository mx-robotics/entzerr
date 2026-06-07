#pragma once
#include <gtk/gtk.h>

void center_offset(int alloc_w, int alloc_h, int img_w, int img_h,
                   double &ox, double &oy);

GdkPixbuf *perspective_warp(GdkPixbuf *src, double pts[4][2]);
