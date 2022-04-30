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


def eval_Bezier1(P, t):
    # P (t) = (1-t) P[0] + tP[1]
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = (1-t) * P[0][xy] + t*P[1][xy]
    return res


def eval_Bezier2(P, t):
    # P(t) = (1-t)^2 * P[0] + 2t(1-t)P[1] + t^2*P[2]
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = (1-t)**2 * P[0][xy] + 2 * t * \
            (1 - t)*P[1][xy] + t**2 * P[2][xy]
    return res


def eval_Bezier3(P, t):
    # P(t) = (1-t)^3 * P[0] + 3(1-t)^2tP[1] + 3(1-t)t^2*P[2] + t^3P[3]
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = (1-t)**3 * P[0][xy] + 3*(1-t)**2 * t * \
            P[1][xy] + 3*(1-t) * t**2 * P[2][xy] + t ** 3 * P[3][xy]
    return res

def eval_dBezier2(P, u, v_factor):
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = ((2*P[0][xy] - 4*P[1][xy] + 2*P[2][xy])
                   * u - 2*P[0][xy] + 2*P[1][xy]) / v_factor
    return res

def eval_dBezier3(P, t, v_factor):
    res = [0.0, 0.0]
    for xy in range(2):
        res[xy] = (-3*(1-t)**2*P[0][xy] + (-6*(1-t)*t + 3*(1-t)**2) * \
            P[1][xy] + (-3*t**2 + 6*(1-t)*t)*P[2][xy] + 3*t**2*P[3][xy]) / v_factor
    return res

def draw_Bezier(P, nsteps):
    xi = P[0][0]
    yi = P[0][1]
    t_delta = 1/nsteps
    t = t_delta
    for ti in range(nsteps):
        if (len(P) == 2):
            p = eval_Bezier1(P, t)
        elif (len(P) == 3):
            p = eval_Bezier2(P, t)
        elif (len(P) == 4):
            p = eval_Bezier3(P, t)

        draw_line(canvas, xi, yi, p[0], p[1], rgb_col(255, 0, 0))
        draw_small_square(xi, yi, rgb_col(255, 255, 0))
        xi = p[0]
        yi = p[1]
        t += t_delta
    for i in range(len(P)):
        draw_small_square(P[i][0], P[i][1], rgb_col(0, 255, 0))

def arc_length_exact(t):
    root = math.sqrt(4*t**2+1)

def do_animation(t):
    global animation_done
    global V_pos
    global V_vec

    v_factor = 5  # reparameterization
    u = t/v_factor
    if (t > v_factor):  # animation stops at t = v_factor
        animation_done = True
    else:
        #B2[1][0] = 7.0 - t
        #B2[1][1] = 0.0 + 3*t/v_factor

        # (2 p_0 - 4 p_1 + 2 p_2) t - 2 p_0 + 2 p_1

        V_vec = eval_dBezier3(B3, u, v_factor)
        V_pos = eval_Bezier3(B3, u)

def f(t):
    p = eval_dBezier3(B3, t, 1)
    return math.sqrt(p[0]**2 + p[1]**2)

def generate_Riemann_reg_grid(x, n, a, b):
    delta_x = (b-a)/n
    for i in range(n+1):
        x.append(a + i*delta_x)

def Simpson_sum(x, a, b):
    s = 0
    n = len(x) - 1
    delta_x_div_6 = (b-a) / (6*n) # assumes equidistant x[i] grid
    for i in range(n):
        s += (f(x[i]) + 4*f((x[i] + x[i+1])/2) + f(x[i+1])) * delta_x_div_6
    return s        

def Simpson_sum_better(int_left, int_right):
    s = 0
    n = n_grid
    delta_x = (int_right - int_left)/n # assumes equidistant x[i] grid
    xi = int_left
    for i in range(n):
        s += (f(xi) + 4*f((xi + xi + delta_x)/2) + f(xi + delta_x))
        xi += delta_x
    return s*delta_x/6

def arc_len_Bezier2_approx(P, nSteps):
    arcLen = 0
    deltaT = 1/nSteps
    t = deltaT

    for i in range(nSteps):
        vector = [0.0,0.0]
        res1 = eval_Bezier3(P, t)
        res2 = eval_Bezier3(P, t + deltaT)

        vector[0] = res2[0] - res1[0]
        vector[1] = res2[1] - res1[1]

        arcLen += math.sqrt(vector[0]**2 + vector[1]**2)

        t += deltaT
    return arcLen   

def draw_scene():
    draw_grid(canvas)
    draw_axis(canvas)
    draw_Bezier(B3, 20)
    GREEN = rgb_col(0, 255, 0)
    draw_dot(canvas, V_pos[0], V_pos[1], GREEN)
    draw_line(canvas, V_pos[0], V_pos[1],
              V_pos[0] + V_vec[0], V_pos[1] + V_vec[1], GREEN)


generate_Riemann_reg_grid(t_grid, n_grid, t_min, t_max)

print(f"Simpson arclen: {Simpson_sum(t_grid, t_min, t_max)}")
print(f"Simpson arclen better: {Simpson_sum_better(t_min, t_max)}")
print(f"Brute force arclen: {arc_len_Bezier2_approx(B3, 1000)}")

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
