import sys
fn, bs = sys.argv[1], int(sys.argv[2])
num, den = [], []
for line in open(fn):
    if line.startswith("#"):
        continue
    t = line.split()
    num.append(float(t[4])); den.append(float(t[5]))
n = len(num)
N, D = sum(num), sum(den)
theta = N / D
blocks = [(sum(num[i:i + bs]), sum(den[i:i + bs]), len(num[i:i + bs])) for i in range(0, n, bs)]
g = len(blocks)
thetas = [(N - bn) / (D - bd) for bn, bd, m in blocks]
theta_jack = g * theta - sum((1 - m / n) * tj for (bn, bd, m), tj in zip(blocks, thetas))
var = 0.0
for (bn, bd, m), tj in zip(blocks, thetas):
    h = n / m
    tau = h * theta - (h - 1) * tj
    var += (tau - theta_jack) ** 2 / (h - 1)
var /= g
print("%g %g" % (theta, var ** 0.5))
