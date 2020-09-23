import sys

imols = [int(x) for x in sys.argv[1:]]
maxmol = max(imols)
nmol = 1
for line in sys.stdin:
    if nmol in imols:
        print (line, end="")
    if line.startswith('$$$$'):
        nmol += 1
        if nmol > maxmol: break