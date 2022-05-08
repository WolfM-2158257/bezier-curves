from tkinter import Tk, Canvas
from graphics_template import *
import math
import time

vp_width, vp_height = 1024, 768
w_xmin, w_ymin, w_xmax = -3, -3, 10
w_ymax = w_ymin + (w_xmax - w_xmin)/vp_width * vp_height
animation_done = False

B1 = [[0.0, 0.0],  # p0
      [7.0, 0.0],  # p1
      [1.0, 4.0]]  # p2

B2 = [[0.0, 0.0],  # p0
      [7.0, 0.0],  # p1
      [1.0, 4.0]]  # p2

B3 = [[0.0, 0.0],  # p0
      [0.5, -0.5],  # p1
      [5.5, -0.5],  # p2
      [7.0, 4.0]]  # p3

V_vec = []
V_pos = []

t_min = 0
t_max = 1
t_grid = []
n_grid = 100

window = Tk()
canvas = Canvas(window, width=vp_width, height=vp_height, bg=rgb_col(0, 0, 0))
canvas.pack()

init_graphics(vp_width, vp_height, w_xmin, w_ymin, w_xmax)


def draw_small_square(xc, yc, fill_col):
    size = 0.05
    draw_line(canvas, xc-size, yc-size, xc+size, yc-size, fill_col)
    draw_line(canvas, xc-size, yc+size, xc+size, yc+size, fill_col)
    draw_line(canvas, xc-size, yc-size, xc-size, yc+size, fill_col)
    draw_line(canvas, xc+size, yc-size, xc+size, yc+size, fill_col)


def init_scene():
    V_pos.append(0.0)
    V_pos.append(0.0)
    V_vec.append(0.0)
    V_vec.append(0.0)
    do_animation(0.0)
    draw_scene()


def eval_Bezier3(P, t):
    # P(t) = (1-t)^3 * P[0] + 3(1-t)^2tP[1] + 3(1-t)t^2*P[2] + t^3P[3]
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = ((1-t)**3 * P[0][xy] + 3*(1-t)**2 * t *
                   P[1][xy] + 3*(1-t) * t**2 * P[2][xy] + t ** 3 * P[3][xy])
    return res


def eval_dBezier3(P, t):
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = (-3*(1-t)**2*P[0][xy] + (-6*(1-t)*t + 3*(1-t)**2) *
                   P[1][xy] + (-3*t**2 + 6*(1-t)*t)*P[2][xy] + 3*t**2*P[3][xy])
    return res


def draw_Bezier(P, nsteps):
    xi = P[0][0]
    yi = P[0][1]
    t_delta = 1/nsteps
    t = t_delta
    for ti in range(nsteps):
        p = eval_Bezier3(P, t)

        draw_line(canvas, xi, yi, p[0], p[1], rgb_col(255, 0, 0))
        draw_small_square(xi, yi, rgb_col(255, 255, 0))
        xi = p[0]
        yi = p[1]
        t += t_delta
    for i in range(len(P)):
        draw_small_square(P[i][0], P[i][1], rgb_col(0, 255, 0))


def d_arc_length_exact(t):
    root = math.sqrt(4*t**2+1)
    return root


def arc_length_exact(t):
    root = math.sqrt(4*t**2+1)
    return 1/2*root*t + 1/4*math.log(root + 2*t)


# de snelheid waarvoor we de integraal berekenen om de lengte van de kromme te berekenen (g'(t) bij get_NR_curve)
def f(t):
    p = eval_dBezier3(B3, t)
    return math.sqrt(p[0]**2 + p[1]**2)

# de totale lengte van de kromme berekenen van links naar rechts (g(t) bij get_NR_curve)
def Simpson_sum(int_left, int_right):
    s = 0
    n = n_grid
    delta_x = (int_right - int_left)/n  # assumes equidistant x[i] grid
    xi = int_left
    for i in range(n):
        s += (f(xi) + 4*f((xi + xi + delta_x)/2) + f(xi + delta_x))
        xi += delta_x
    return s*delta_x/6


# gebruiken newton raphson om g(t) - s = 0 te benaderen ()
def get_NR_curve_t(s):
    # set t0
    t = t_min + s/Simpson_sum(t_min, t_max)*(t_max-t_min)
    # print("init t0", t)
    for it in range(5):
        t -= (Simpson_sum(t_min, t)-s)/f(t)
        # print("it t", t)
    return t


def pythagoras(p):
    res = math.sqrt(p[0]**2+p[1]**2)
    return res


def do_animation(s):
    global animation_done
    global V_pos
    global V_vec

    speed = 2
    u = s*speed
    t = get_NR_curve_t(u)  # berekent genormalizeerde t
    if (t > t_max):
        animation_done = True
    else:
        V_pos = eval_Bezier3(B3, t)  # positie

        v = eval_dBezier3(B3, t)  # snelheids vector
        V_vec[0], V_vec[1] = v[0] / pythagoras(v) * speed, v[1]/pythagoras(v) * speed


def generate_Riemann_reg_grid(x, n, a, b):
    delta_x = (b-a)/n
    for i in range(n+1):
        x.append(a + i*delta_x)


def draw_scene():
    draw_grid(canvas)
    draw_axis(canvas)
    draw_Bezier(B3, 20)
    GREEN = rgb_col(0, 255, 0)
    draw_dot(canvas, V_pos[0], V_pos[1], GREEN)
    draw_line(canvas, V_pos[0], V_pos[1],
              V_pos[0] + V_vec[0], V_pos[1] + V_vec[1], GREEN)


init_time = time.perf_counter()
prev_draw_time = 0
init_scene()
while (not animation_done):
    draw_dt = time.perf_counter() - init_time - prev_draw_time
    if (draw_dt > 0.02):  # 50 fps
        prev_draw_time += draw_dt
        do_animation(prev_draw_time)
        canvas.delete("all")
        draw_scene()
        canvas.update()
