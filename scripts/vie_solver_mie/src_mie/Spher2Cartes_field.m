function [StoC] = Spher2Cartes_field(theta,phi)


StoC = [ sin(theta).*cos(phi)   cos(theta).*cos(phi) -sin(phi)
         sin(theta).*sin(phi)   cos(theta).*sin(phi)  cos(phi)
         cos(theta)           -sin(theta)           0*theta     ];