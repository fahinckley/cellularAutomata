%--------------------------------
% Draws a NACA standard series airfoil for the specified parameter values
%--------------------------------
% Franklin Hinckley
% 11 April 2016
%--------------------------------
% Inputs:
%   p1 (double):
%   p2 (double):
%   p3 (double):
%   c (double): chord length [m]
%   numSeg (double): number of segments
% Outputs:
%   xs (numSeg x 1 double array): x-coordinates of segment endpoints
%   ys (numSeg x 1 double array): y-coordinates of segment endpoints
%--------------------------------

function [xs,ys] = NACA(p1,p2,p3,AoA,c,numSeg)

% Compute intermediate values
m = p1/100; % maximum camber
p = p2/10;  % location of maximum camber
t = p3/100; % thickness as a fraction of chord
AoA = -AoA; % fix convention for angle of attack

% Determine length array using cosine spacing
beta = linspace(0,pi,numSeg);
x = ((1-cos(beta))./2)*c;

yt = (t/.2)*c.*((.2969*sqrt(x./c))-(.1260*(x./c))-(.3516*((x./c).^2))+(.2843*((x./c).^3))-...
    (.1036*((x./c).^4)));

xc = x*cosd(AoA);
yc = x*sind(AoA);

% Determine if symmetric or not and determine shape as appropriate
if p1==0 && p2==0
    xu = xc - (yt*sind(AoA));
    xl = xc + (yt*sind(AoA));
    yu = yc + (yt*cosd(AoA));
    yl = yc - (yt*cosd(AoA));
else
    xu = zeros(size(x));
    xl = zeros(size(x));
    yu = zeros(size(x));
    yl = zeros(size(x));
    dycdx = zeros(size(x));
    thetau = zeros(size(x));
    thetal = zeros(size(x));
    for ii=1:length(x)
        if x(ii)<=(p*c)
            yc(ii)=m*(x(ii)/(p^2))*((2*p)-(x(ii)/c));
            %yc=(m/(p^2))*(2*p.*x-(x.^2));
            dycdx(ii)=((2*m)/(p^2))*(p-(x(ii)));
            thetau(ii)=atan(dycdx(ii));
            thetal(ii)=-atan(dycdx(ii));
        else
            yc(ii)=m*((c-x(ii))/((1-p)^2))*(1+(x(ii)/c)-(2*p));
            %yc=(m/((1-p)^2))*(1-(2*p)+(2*p.*x)-(x.^2));
            dycdx(ii)=((2*m)/((1-p)^2))*(p-(x(ii)));
            thetau(ii)=atan(dycdx(ii));
            thetal(ii)=-atan(dycdx(ii));
        end
    xu(ii) = x(ii)  - (yt(ii)*sin(thetau(ii)));
    xl(ii) = x(ii)  + (yt(ii)*sin(thetal(ii)));
    yu(ii) = yc(ii) + (yt(ii)*cos(thetau(ii)));
    yl(ii) = yc(ii) - (yt(ii)*cos(thetal(ii)));
    end
end

%Combine the upper and lower surface series
xsl=fliplr(xl(1:end-1));
xs=[xu xsl];
ysl=fliplr(yl(1:end-1));
ys=[yu ysl];

xs(end)=xs(1);
ys(end)=ys(1);

xscl=fliplr(xc(1:end-1));
xsc=[xc xscl];
yscl=fliplr(yc(1:end-1));
ysc=[yc yscl];

%Plot the airfoil
hold on
plot(xs,ys)
plot(xc,yc,'r')
axis equal
hold off
drawnow

end