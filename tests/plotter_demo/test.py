import taichi as ti
ti.init(print_ir = True)
n = 10
a = ti.field(ti.f32, shape = (n))

@ti.kernel
def fun(m: ti.i32):
    a[9] = m
    for i in range(m):
        a[i] = i

fun(7)
print(a.to_numpy().sum())
