function ipi=ipf(cr,D,w)
% ipf    Calculate ideal particle image. 
% Usage: ipi=ipf(cr,D,w)
%
% Calculates an ideal particle image ipi.  The particle has diameter D and
% width parameter w.  2w is the width of 76% of the fall off.

% revision history:
% 08/04/05 Mark D. Shattuck <mds> ipf.m  
% 01/30/06 mds added abs(cr)
% 04/30/07 mds made w a true measure of width.

ipi=(1-tanh((abs(cr)-D/2)/w))/2;