function flag = isInsideContour(z, center, radius, varargin)
% This short scripts checks wheter the points z are inside the circle of
% center "center" and radius "radius". The varargin is used in the
% elliptical case, where "str" is the ratio vertical axis/horizontal axis
% and "rot" is the clockwise rotation from the Y-axis.

str = 1;
rot = 0;
if nargin >= 4
    str = varargin{1};
    if nargin >= 5
        rot = varargin{2};
    end
end
  w = (z - center)/radius;
  flag = (real(w)*cos(rot) + imag(w)*sin(rot)).^2 + (real(w)*sin(rot) - imag(w)*cos(rot)).^2/str^2 <= 1;
  end
