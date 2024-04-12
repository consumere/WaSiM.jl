#Mandelbrot test
using PyPlot

function mandelbrot(width, height, max_iter)
    img = zeros(UInt8, height, width, 3)
    x_range = range(-2, stop=1, length=width)
    y_range = range(-1.5, stop=1.5, length=height)
    
    for j = 1:height
        for i = 1:width
            c = x_range[i] + y_range[j]*im
            z = 0.0
            iter = 0
            
            while abs(z) <= 2 && iter < max_iter
                z = z^2 + c
                iter += 1
            end
            
            if iter == max_iter
                img[j, i, :] = [0, 0, 0]  # Black for points in the Mandelbrot set
            else
                img[j, i, :] = [iter % 256, iter % 256, iter % 256]  # Color based on iteration count
            end
        end
    end
    
    return img
end

function plot_mandelbrot(image)
    height, width, _ = size(image)
    x_range = range(-2, stop=1, length=width)
    y_range = range(-1.5, stop=1.5, length=height)
    
    fig = figure(figsize=(8, 6))
    ax = fig[:add_subplot](111)
    ax[:set_aspect]("equal")
    ax[:set_xlim](-2, 1)
    ax[:set_ylim](-1.5, 1.5)
    ax[:axis]("off")
    
    ax[:imshow](image[:, :, 1], extent=[-2, 1, -1.5, 1.5], cmap="viridis", alpha=0.7)
    
    show()
end

ig = mandelbrot(800, 600, 100);
#ig = mandelbrot(33, 99, 3333);
plot_mandelbrot(ig)