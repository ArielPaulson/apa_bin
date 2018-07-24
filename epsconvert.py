import Image
import sys
s = sys.argv[1] 
v = sys.argv[2] + ".png"
im = Image.open(s)
im.save(v) 

