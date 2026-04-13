import matplotlib.pyplot as plt
import control as ct

# (a)
sys_a = ct.tf([1], [1, 4, 3])
print("System (a):")
print(sys_a)

t, y = ct.step_response(sys_a)
plt.figure()
plt.plot(t, y)
plt.title('(a) Step Response')
plt.grid()

# (b)
sys_b = ct.tf([1], [1, 4, 4])
print("System (b):")
print(sys_b)

t, y = ct.step_response(sys_b)
plt.figure()
plt.plot(t, y)
plt.title('(b) Step Response')
plt.grid()

# (c)
sys_c = ct.tf([1], [1, 4, 8])
print("System (c):")
print(sys_c)

t, y = ct.step_response(sys_c)
plt.figure()
plt.plot(t, y)
plt.title('(c) Step Response')
plt.grid()

# (d)
sys_d = ct.tf([1], [1, 0, 36])
print("System (d):")
print(sys_d)

t, y = ct.step_response(sys_d)
plt.figure()
plt.plot(t, y)
plt.title('(d) Step Response')
plt.grid()

plt.show()