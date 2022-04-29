# graphics template Wisk I
import math

WIDTH, HEIGHT = 1024, 768

window_xmin = -1.0
window_ymin = -3.0
window_xmax = 7.0
window_ymax = window_ymin + (window_xmax - window_xmin)/WIDTH * HEIGHT
# above calculation to give window the same aspect ratio as viewport

def init_graphics (w, h, w_xmin, w_ymin, w_xmax):
    global WIDTH, HEIGHT, window_xmin, window_ymin, window_xmax, window_ymax
    WIDTH = w
    HEIGHT = h
    window_xmin = w_xmin
    window_ymin = w_ymin
    window_xmax = w_xmax
    window_ymax = window_ymin + (window_xmax - window_xmin)/WIDTH * HEIGHT
    # above calculation to give window the same aspect ratio as viewport


def rgb_col (r, g, b):
    """
    Assumes r,g and b integer values 0..255
    Returns hex_coded color string '#rrggbb'
    """
    c = '0123456789abcdef'
    return '#' + c[int(r/16)] + c[r%16] + c[int(g/16)] + c[g%16] + c[int(b/16)] + c[b%16]


def window_to_viewport (p_w):
    # translate point over (-window_xmin, -window_ymin)
    x_vp = p_w[0] - window_xmin
    y_vp = p_w[1] - window_ymin

    # scale point with
    # (viewport width / window width , viewport height / window height)
    x_vp = x_vp / (window_xmax - window_xmin) * WIDTH
    y_vp = y_vp / (window_ymax - window_ymin) * HEIGHT

    # flip vertical axis and translate to put (0,0) at lower left corner
    y_vp = HEIGHT - y_vp
    
    return [x_vp,y_vp]


def draw_pixel (canvas, xp, yp, fill_col):
    vp_point = window_to_viewport([xp,yp])
    canvas.create_line(vp_point[0], vp_point[1], vp_point[0]+1, vp_point[1]+1,
                       fill = fill_col)


def draw_dot (canvas, xd, yd, fill_col):
    vp_dot = window_to_viewport([xd,yd])
    r = 5 # radius of the dot in pixels
    canvas.create_oval(vp_dot[0]-r, vp_dot[1]-r, vp_dot[0]+r, vp_dot[1]+r,
                       fill=fill_col)


def draw_line (canvas, x0, y0, x1, y1, fill_col):
    vp_p0 = window_to_viewport([x0,y0])
    vp_p1 = window_to_viewport([x1,y1])
    canvas.create_line(vp_p0[0], vp_p0[1], vp_p1[0], vp_p1[1],
                       fill = fill_col)


def draw_rect (canvas, x0, y0, x1, y1, fill_col, to_be_filled):
    vp_p0 = window_to_viewport([x0,y0])
    vp_p1 = window_to_viewport([x1,y1])
    if to_be_filled:
        canvas.create_rectangle(vp_p0[0], vp_p0[1], vp_p1[0], vp_p1[1],
                                fill=fill_col, outline=fill_col)
    else:
        canvas.create_rectangle(vp_p0[0], vp_p0[1], vp_p1[0], vp_p1[1],
                                outline=fill_col)
        

def draw_circle (canvas, xc, yc, r, fill_col, to_be_filled):
    vp_c = window_to_viewport([xc,yc])
    vp_r = window_to_viewport([0,r])
    vp_0 = window_to_viewport([0,0])
    r = vp_r[1] - vp_0[1]
    if to_be_filled:
        canvas.create_oval(vp_c[0]-r, vp_c[1]-r, vp_c[0]+r, vp_c[1]+r,
                           fill= fill_col)
    else:
        canvas.create_oval(vp_c[0]-r, vp_c[1]-r, vp_c[0]+r, vp_c[1]+r,
                           outline=fill_col)


def draw_axis (canvas):
    draw_line (canvas, 0, window_ymin, 0, window_ymax, rgb_col(255,255,255))
    draw_line (canvas, window_xmin, 0, window_xmax, 0, rgb_col(255,255,255))
    

def draw_grid (canvas):
    for i in range (math.trunc(window_xmin), math.ceil(window_xmax)):
        draw_line (canvas, i, window_ymin, i, window_ymax, rgb_col(0,0,255))
    for j in range (math.trunc(window_ymin), math.ceil(window_ymax)):
        draw_line (canvas, window_xmin, j, window_xmax, j, rgb_col(0,0,255))

