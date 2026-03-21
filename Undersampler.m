classdef Undersampler
   properties    
      
   end
   methods
       function obj = Undersampler

       end

     
       function [Ku, U] = undersample(us, K, accel, method)
           [Nx,Ny,Nt] = size(K);
           U = false(size(K));
           switch method
               case 'poisson_disc'
                    for nt = 1:Nt
                        U(:,:,nt) = poisson_disc(Nx,Ny,accel,20,0.05);
                    end
           end
           Ku = U .* K;
           % figure();
           % nslices = 3;
           % indices = round(linspace(1,size(U,3),nslices));
           % slice(double(permute(U,[2 3 1])), indices,[], []);
           drawnow;
       end

   end

end

function u = poisson_disc(nx, ny, accel, sp, tol)
    x = -(nx-1)/2:(nx-1)/2;
    y = -(ny-1)/2:(ny-1)/2;
    [x,y] = ndgrid(x,y);
    r = sqrt(x.^2 + y.^2);
    smin = 0;
    smax = max(nx,ny);
    while smin < smax
        s = (smax + smin) / 2;
        Rx = (1 + s*r) * nx / max(nx, ny);
        Ry = (1 + s*r) * ny / max(nx, ny);
        u = sample(Rx, Ry, sp, nx, ny);
        actual_accel = (nx * ny) / nnz(u);
        if abs(actual_accel - accel) < tol
            break;
        end
        if actual_accel < accel
            smin = s;
        else
            smax = s;
        end
    end
end

function u = sample(Rx, Ry, sp, nx, ny)
    c = 2;
    lx = ceil((nx+1)/2-c):floor((nx+1)/2+c);
    ly = ceil((ny+1)/2-c):floor((ny+1)/2+c);
    u = false(nx,ny);
    u(lx,ly) = 1;
    % [x, y] = ndgrid(1:nx, 1:ny);
    % cx = (nx+1)/2;
    % cy = (ny+1)/2;
    % u = ((x - cx).^2 + (y - cy).^2) <= c^2;
    pxs = zeros(nx*ny, 1);
    pys = zeros(nx*ny, 1);
    pxs(1) = randi(nx);
    pys(1) = randi(ny);
    num_actives = 1;
    while num_actives >= 1
        i = randi(num_actives);
        px = pxs(i);
        py = pys(i);
        rx = Rx(px,py);
        ry = Ry(px,py);
        done = false;
        k = 1;
        while ~done && k <= sp
            v = sqrt(rand(1) * 3 + 1);
            t = 2 * pi * rand(1);
            qx = px + v * rx * cos(t);
            qy = py + v * ry * sin(t);
            if qx>=1 && qx < nx + 1 && qy >=1 && qy < ny + 1
                startx = max(fix(qx - rx), 1);
                endx = min(fix(qx+rx+1),nx);
                starty = max(fix(qy - ry), 1);
                endy = min(fix(qy+ry+1),ny);
                done = true;
                for x = startx:endx
                    for y = starty:endy
                        if (u(x,y) == 1 && (((qx-x)/Rx(x,y))^2 + ((qy-y)/Ry(x,y))^2 < 1))                            
                            done = false;
                            break;
                        end
                    end
                end
            end
            k = k+1;
        end
        if done
            num_actives = num_actives + 1;
            pxs(num_actives) = fix(qx);
            pys(num_actives) = fix(qy);
            u(fix(qx), fix(qy)) = 1;
        else
            pxs(i) = pxs(num_actives);
            pys(i) = pys(num_actives);
            num_actives = num_actives - 1;
        end
    end
end