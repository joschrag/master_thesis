function c = valid_eq(p,mode,opt)
arguments
    p {mustBePrimeOrZero} = 0;
    mode {mustBeTextScalar} = "3Q3"
    opt.verbose {mustBeInRange(opt.verbose,0,2)} = 1;
    opt.stream = rng("default");
end
if p==0
    x0 = 5*rand(opt.stream)-2.5;
    y0 = 5*rand(opt.stream)-2.5;
    z0 = 5*rand(opt.stream)-2.5;
else
    x0 = randi(p-1);
    y0 = randi(p-1);
    z0 = randi(p-1);
end
if opt.verbose > 0
fprintf("Generating %s problem with intersection in (%f,%f,%f).\n",mode,x0,y0,z0)
end
if strcmp(mode,"3Q3")
    x_idx = [1,7,10];
    y_idx = [4,8,10];
    z_idx = [6,9,10];
    xy_idx = [2,7,8,10];
    xz_idx = [3,7,9,10];
    yz_idx = [5,8,9,10];
    xv = [1,-x0];
    yv = [1,-y0];
    zv = [1,-z0];
    c = zeros(3,10);
    c(:,[7,10]) =c(:,[7,10])+ [1,-x0];
    c(:,[8,10]) = c(:,[8,10])+[1,-y0];
    c(:,[9,10]) = c(:,[9,10])+[1,-z0];
    c(:,x_idx) = c(:,x_idx) + [1;4;2].*conv(xv,xv);
    c(:,y_idx) = c(:,y_idx) + [5;4;3].*conv(yv,yv);
    c(:,z_idx) = c(:,z_idx) + [1;2;1].*conv(zv,zv);
    c(:,xy_idx) = c(:,xy_idx) + [2;2;5].*reshape((yv.*xv.')',1,[]);
    c(:,xz_idx) = c(:,xz_idx) + [3;5;7].*reshape((zv.*xv.')',1,[]);
    c(:,yz_idx) = c(:,yz_idx) + [4;-3;-1].*reshape((zv.*yv.')',1,[]);
elseif strcmp(mode,"3C3")
    x_idx = [1,11,17,20];
    y_idx = [7,14,18,20];
    z_idx = [10,16,19,20];
    xxy_idx = [2,11,12,17,18,20];
    xxz_idx = [3,11,13,17,19,20];
    xyy_idx = [4,14,12,18,17,20];
    xzz_idx = [6,16,13,19,16,20];
    yyz_idx = [8,14,15,18,19,20];
    yzz_idx = [9,16,15,19,18,20];
    xyz_idx = [5,12,11,17,15,19,18,20];
    xy_idx = [12,17,18,20];
    xz_idx = [13,17,19,20];
    yz_idx = [15,18,19,20];
    xv = [1,-x0];
    yv = [1,-y0];
    zv = [1,-z0];
    c = zeros(3,20);
    c(:,x_idx) = c(:,x_idx) + [1;7;2].*conv(conv(xv,xv),xv);
    c(:,y_idx) = c(:,y_idx) + [5;1;-2].*conv(conv(yv,yv),yv);
    % c(:,z_idx) = c(:,z_idx) + [1;-4;-2].*conv(conv(zv,zv),zv);
    % c(:,xxy_idx) = c(:,xxy_idx) + [2;4;-2].*reshape((yv.*conv(xv,xv).')',1,[]);
    % c(:,xxz_idx) = c(:,xxz_idx) + [-2;1;3].*reshape((zv.*conv(xv,xv).')',1,[]);
    % c(:,xyy_idx) = c(:,xyy_idx) + [7;-3;-2].*reshape((xv.*conv(yv,yv).')',1,[]);
    % c(:,xzz_idx) = c(:,xzz_idx) + [11;-3;2].*reshape((xv.*conv(zv,zv).')',1,[]);
    c(:,yyz_idx) = c(:,yyz_idx) + [-11;1;5].*reshape((zv.*conv(yv,yv).')',1,[]);
    % c(:,yzz_idx) = c(:,yzz_idx) + reshape((yv.*conv(zv,zv).')',1,[]);
    % c(:,xyz_idx) = c(:,xyz_idx) + reshape((zv.*reshape((yv.*xv.')',1,[]).')',1,[]);
    xy = reshape((yv.*xv.')',1,[]);
    c(:,xy_idx) = c(:,xy_idx) -(c(:,xy_idx(1))./xy(1)).*xy ;
    yz = reshape((zv.*yv.')',1,[]);
    c(:,yz_idx) = c(:,yz_idx) + [2;-3;-2].*yz;
    xz = reshape((zv.*xv.')',1,[]);
    c(:,xz_idx) = c(:,xz_idx) -(c(:,xz_idx(1))./xz(1)).*xz ;
    c(:,x_idx(end-2:end)) = c(:,x_idx(end-2:end))+conv(xv,xv);
    y_quad = conv(yv,yv);
    c(:,y_idx(end-2:end)) = c(:,y_idx(end-2:end))-(c(:,y_idx(end-2))./y_quad(1)).*y_quad;
    z_quad = conv(zv,zv);
    c(:,z_idx(end-2:end)) = c(:,z_idx(end-2:end))-(c(:,z_idx(end-2))./z_quad(1)).*z_quad;
    c(:,x_idx(end-1:end)) = c(:,x_idx(end-1:end))+xv;
    c(:,y_idx(end-1:end)) = c(:,y_idx(end-1:end))-(c(:,y_idx(end-1))./yv(1)).*yv;
    c(:,z_idx(end-1:end)) = c(:,z_idx(end-1:end))-(c(:,z_idx(end-1))./zv(1)).*zv;
end
c = mod(c,p);
control = mod(c*mod([x0^3;x0^2*y0;x0^2*z0;x0*y0^2;x0*y0*z0;x0*z0^2;y0^3;y0^2*z0;y0*z0^2;z0^3;x0^2;x0*y0;x0*z0;y0;y0*z0;z0^2;x0;y0;z0;1],p),p);
if opt.verbose > 0
    fprintf("Control: (%f,%f,%f).\n",control)
end
end