import re

def switch_class(classes,to_add,all_classes):
    classes = set(classes.split(" "))
    for c in all_classes:
        if c in classes:
            classes.remove(c)
    for c in to_add:
        classes.add(c)
    
    return " ".join(classes)
#https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color

def get_rgb(color_string:str):
    if color_string.startswith("#"):
        r=int(color_string[1:3],16)/255
        g=int(color_string[3:5],16)/255
        b=int(color_string[5:7],16)/255
    elif color_string.startswith("rgb"):
        m = re.match(r'rgb\(([0-9]+), ?([0-9]+), ?([0-9]+)\)',color_string)
        r,g,b = list(map(lambda x :int(x)/255,m.groups()))
    return r,g,b
def get_L(color_string):
    r,g,b = get_rgb(color_string)
    L = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return L
#https://github.com/Myndex/apca-w3/?tab=readme-ov-file
def get_Lc(S_apc):
    return S_apc*100
def get_S_apc(C,P_out=0.1,W_offset=0.027):
    if(abs(C)<P_out):
        return 0
    if C>0:
        return C-W_offset
    else:
        return C+W_offset

def get_C(y_bg,Y_txt,R_scale=1.14,P_in=0.0005):
    if(abs(y_bg-Y_txt)<P_in):
        return 0
    if Y_txt<y_bg:
        return get_S_norm(y_bg,Y_txt)*R_scale
    return get_S_rev(y_bg,Y_txt)*R_scale

def get_S_norm(y_bg,Y_txt):
    return y_bg**0.56-Y_txt**0.57
    
def get_S_rev(y_bg,Y_txt):
    return y_bg**0.65-Y_txt**0.62
    
def Ys(r,g,b):
    return f_clamp(((r)**2.4)*0.2126729+((g)**2.4)*0.7151522+((b)**2.4)*0.0721750)
def f_clamp(y,b_thr=0.022,b_exp=1.414):
    print(y)
    if y>b_thr:
        return y
    return y+(b_thr-y)*b_exp

def get_text_color(background_color):
    y_bg = Ys(*get_rgb(background_color))
    y_white = Ys(1,1,1)
    y_black = Ys(0,0,0)
    l_white = get_Lc(get_S_apc(get_C(y_bg,y_white)))
    l_black = get_Lc(get_S_apc(get_C(y_bg,y_black)))
    if(abs(l_white)>abs(l_black)):
        return "#ffffff"
    return "#000000"
    # L = get_L(background_color)
    # print(background_color,L)
    # if L < 0.4:
    #     return "#ffffff"
    # return "#000000"