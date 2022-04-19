function F = FIcurve(I, varargin)
    %F(I) for vector I; a gain factor, d is noise factor; 
    p = inputParser;
    p.addParameter('a', 270);
    p.addParameter('b', 108);
    p.addParameter('d', 0.154);
    p.parse(varargin{:});
    Itemp = p.Results.a*I - p.Results.b;
    F = Itemp./(1.-exp(-p.Results.d*(Itemp)));
    F(Itemp == 0) = 1/p.Results.d;
end

