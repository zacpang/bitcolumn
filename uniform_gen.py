import random
f = open("7_differ_attr_new.txt","w");
start = 1;
enda = pow(2,32)-1;
endc = pow(2,22);
ende = pow(2,25);
endn = pow(2,19);
for i in range(1,20000000):
        b = 0;d = 0;m = 0;
        a =int(random.uniform(start, enda));
        while(b <= 0):
                b =int(random.gauss(pow(2,25),pow(2,25)));
        c =int(random.uniform(start, endc));
        while(d <= 0):
                d =int(random.gauss(pow(2,23),pow(2,23)));
        e =int(random.uniform(start, ende));
        while(m <= 0):
                m =int(random.gauss(pow(2,26),pow(2,26)));
        n =int(random.uniform(start, endn));
        f.write(str(a)+','+str(b)+','+str(c)+','+str(d)+','+str(e)+','+str(m)+','+str(n));
        f.write('\n');
f.close();
