def Euler(init_x, y, step, x, function):

    while init_x < x:
        y = y + step*function(x,y) 
        init_x = init_x + step

    return y