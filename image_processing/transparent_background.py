import sys
from PIL import Image


if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <image file>' %sys.argv[0]
	exit()

imgfile = sys.argv[1]
outimgfile = imgfile[:imgfile.find('.')]+'_trans.png'

img = Image.open(imgfile)
img = img.convert("RGBA")
datas = img.getdata()


# grey pixel colors 
grey = [224,192,160,128]

newData = []
for item in datas:
	# # change all pixels that dont have an pixel value of 0 to transparent
	# if item[0] or item[1] or item[2]:
	# 	newData.append((255, 255, 255, 0))

	# change white pixels to transparent
	if item[0] == 255 and item[1] == 255 and item[2] == 255:
		newData.append((255, 255, 255, 0))

	else:
		newData.append(item)

img.putdata(newData)
img.save(outimgfile, "PNG")
