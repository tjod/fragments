import sys

imol = int(sys.argv[1])
nmol = 1
for line in sys.stdin:
    if nmol == imol:
        print (line, end="")
    if line.startswith('$$$$'):
        nmol += 1
        if nmol > imol: break