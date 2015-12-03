function Intens = BC(lambda,theta);% BC   Calculates the intensity distribution%    for a blazed reflection grating of up to 50 slits.  % -----------------------------------------------% Set up the optical and mathematical variables  % -----------------------------------------------% All optical variables are global and set elsewhere  global blazAngl stepLeng numSlits;  % -------------------------------------------------% Calculate the grating intensity% -------------------------------------------------% First calculate alpha and beta.% Note: the vector theta is calculated by the Driver.  phi  = pi * stepLeng * sin(theta) / lambda;% Note: beta is calculated as in a transmission grating% except that (1) the effective slit width is equal to % the effective "slit separation, and (2) the "single% slit" pattern is centered around the angle 2*blazAngl  beta = pi * stepLeng .* sin(theta-pi*blazAngl/90) / lambda;% First calculate the diffraction pattern from a single slit  slitPatn = (sin(beta+eps)./(beta+eps)).^2;  % Then calculate the diffraction pattern from N point sources  srcePatn = (sin(numSlits*(phi+eps))./(numSlits*sin(phi+eps))).^2;% Return the two intensity patterns separately  Intens = [slitPatn;srcePatn];  return;