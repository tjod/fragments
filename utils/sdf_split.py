import sys

if len(sys.argv) < 5:
    print ("usage: sdf_split input output1 output2 mod_n")
    sys.exit()
infile = sys.argv[1]
out1 = sys.argv[2]
out2 = sys.argv[3]
mod_n = int(sys.argv[4])
fp1 = open(out1, "w")
fp2 = open(out2, "w")
m = []
nmol = 0
with open(infile) as fp:
    for line in fp:
        m.append(line)
        if line.startswith('$$$$'):
            if nmol % mod_n == 0:
                fp2.write("".join(m))
            else:
                fp1.write("".join(m))
            m = []
            nmol += 1