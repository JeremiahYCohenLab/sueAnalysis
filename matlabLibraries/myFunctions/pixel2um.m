function um = pixel2um(pixel, mag)
%convert pixel length to length in um. calculated by Sue 12/09/2020
um = 1000 * pixel/1.3132 * 4/mag;