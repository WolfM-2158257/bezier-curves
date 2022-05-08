from tkinter import Tk, Canvas
from graphics_template import *
import math, time

vp_width, vp_height = 1024, 768
w_xmin, w_ymin, w_xmax = -3, -3, 10
w_ymax = w_ymin + (w_xmax - w_xmin)/vp_width * vp_height
animation_done = False

B1 =    [[0.0, 0.0], # p0
        [7.0, 0.0],  # p1
        [1.0, 4.0]]  # p2

window = Tk()
canvas = Canvas(window, width = vp_width, height=vp_height, bg=rgb_col(0,0,0))
canvas.pack()

init_graphics(vp_width, vp_height, w_xmin, w_ymin, w_xmax);


def draw_small_square(xc, yc, fill_col):
    size = 0.05
    draw_line(canvas, xc-size, yc-size, xc+size, yc-size, fill_col)
    draw_line(canvas, xc-size, yc+size, xc+size, yc+size, fill_col)
    draw_line(canvas, xc-size, yc-size, xc-size, yc+size, fill_col)
    draw_line(canvas, xc+size, yc-size, xc+size, yc+size, fill_col)


def do_animation (t):
    global animation_done
    duration = 5
    if (t > duration): # animation stops at t = duration
        animation_done=True
    else:
        # transleer 2de punt
        B1[1][0] = 7.0 - t
        B1[1][1] = 0.0 + 6*t/duration



def init_scene():
    # no data inits needed
    do_animation (0.0)
    draw_scene();


def eval_Bezier2(P, t):
    # P(t) = (1-t)^2 * P[0] + 2t(1-t)P[1] + t^2*P[2]
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = (1-t)**2 * P[0][xy] + 2 * t * (1 - t)*P[1][xy] + t**2 * P[2][xy]
    return res

def draw_Bezier(P, nsteps):
    xi = P[0][0]
    yi = P[0] [1]
    t_delta = 1/nsteps
    t = t_delta
    for ti in range (nsteps):
        p = eval_Bezier2(P, t)
        draw_line (canvas, xi, yi, p[0], p[1], rgb_col (255,0,0))
        draw_small_square (xi, yi, rgb_col (255,255,0))
        xi = p [0]
        yi = p [1]
        t += t_delta
    for i in range(len(P)):
        draw_small_square(P[i][0], P[i][1], rgb_col(0,255, 0))


def draw_scene ():
    draw_grid(canvas)
    draw_axis(canvas)
    draw_Bezier(B1, 20)


init_time = time.perf_counter()
prev_draw_time = 0
init_scene()
while (not animation_done):
    draw_dt = time.perf_counter() - init_time - prev_draw_time
    if (draw_dt > 0.02): #50 fps
        prev_draw_time += draw_dt
        do_animation (prev_draw_time)
        canvas.delete("all")
        draw_scene()
        canvas.update()
