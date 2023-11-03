import numpy as np
import cv2 as cv

drawing = False

interv_s = 0.2

x_vals = []
y_vals = []
z_vals = []

# mouse callback function
def draw_circle(event,x,y,flags,param):
    global drawing
    global vals
    if event == cv.EVENT_LBUTTONDOWN:
        drawing = True
        cv.circle(img,(x,y),1,(0,255,0),-1)
        vals.append((-y+256)/16)
    elif event == cv.EVENT_LBUTTONUP:
        drawing = False
        cv.circle(img,(x,y),1,(0,255,0),-1)
        vals.append((-y+256)/16)
    elif drawing and event == cv.EVENT_MOUSEMOVE:
        cv.circle(img,(x,y),1,(0,255,0),-1)
        vals.append((-y+256)/16)

# Create a black image, a window and bind the function to window
img = np.zeros((512,1536,3), np.uint8)
img[255, :, 0] = 0
img[255, :, 1] = 0
img[255, :, 2] = 255
cv.namedWindow('canvas')
cv.setMouseCallback('canvas',draw_circle)
vals = x_vals
while(1):
    cv.imshow('canvas',img)
    if cv.waitKey(20) & 0xFF == 27:
        break
cv.destroyAllWindows()

# Create a black image, a window and bind the function to window
img = np.zeros((512,1536,3), np.uint8)
img[255, :, 0] = 0
img[255, :, 1] = 0
img[255, :, 2] = 255
cv.namedWindow('canvas')
cv.setMouseCallback('canvas',draw_circle)
vals = y_vals
while(1):
    cv.imshow('canvas',img)
    if cv.waitKey(20) & 0xFF == 27:
        break
cv.destroyAllWindows()

# Create a black image, a window and bind the function to window
img = np.zeros((512,1536,3), np.uint8)
img[255, :, 0] = 0
img[255, :, 1] = 0
img[255, :, 2] = 255
cv.namedWindow('canvas')
cv.setMouseCallback('canvas',draw_circle)
vals = z_vals
while(1):
    cv.imshow('canvas',img)
    if cv.waitKey(20) & 0xFF == 27:
        break
cv.destroyAllWindows()

l = min(len(x_vals), len(y_vals), len(z_vals))

f = open('traj_delta4.txt', 'w')

for i in range(l):
    f.write(str(i*0.2) + " " + str(x_vals[i]) + " " + str(y_vals[i]) + " " + str(z_vals[i]) + " 1.0 0.0 0.0 0.0\n")

f.close()
