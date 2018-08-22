fig, axes = plt.subplots(len(periods), 1, figsize=(15, 10))
d = 0.05
dt = trace.stats.delta
ug = trace.data
kg = len(trace.data)
dt_in = trace.stats.delta
for period, sa, ax in zip(periods, sas, axes):
    w = 2 * np.pi / period
    pr = 2 * np.pi/w
    nn=1
    wd=np.sqrt(1.-d*d)*w
    w2=w*w
    w3=w2*w
    f2=1./w2
    f3=d*w
    f4=1./wd
    f5=f3*f4
    f6=2.*f3
    a = np.zeros(kg)
    v = np.zeros(kg)
    dd = np.zeros(kg)

    ns= int(10.*dt_in/pr-0.01)+1
    dt=dt_in/float(ns)

    for k in range(kg):
        f1=2.*d/w3/dt
        e=np.exp(-f3*dt)
        g1=e*dsin(wd*dt)
        g2=e*dcos(wd*dt)
        h1=wd*g2-f3*g1
        h2=wd*g1+f3*g2
        dug=(ug(k+1)-ug(k))/float(ns)
        g=ug(k)
        z1=f2*dug
        z3=f1*dug
        z4=z1/dt

        for i in range(1,ns):
            z2=f2*g
            b=d[i-1]+z2-z3
            a=f4*v[i-1]+f5*b+f4*z4
            dd[i]=a*g1+b*g2+z3-z2-z1
            v[i]=a*h1-b*h2-z4
            a[i]=-f6*v[i]-w2*d[i]
            nn = nn + 1
            g=g+dug
    ax.plot(trace.times(), a, label='SA'  +  str(period)  +  ' [cm / s / s]')
    ax.plot(trace.times(), trace.data, label='Acceleration [cm / s / s]')
    maxy = np.max(np.abs(a))
    print(maxy, sa, np.abs(maxy  -  sa) /  ((maxy  +  sa)  *  0.5)  *  100)
    legend = ax.legend()
