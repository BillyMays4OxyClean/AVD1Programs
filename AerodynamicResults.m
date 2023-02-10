classdef AerodynamicResults
    properties
        AoA {mustBeNumeric}
        CL {mustBeNumeric}
        CD {mustBeNumeric}
        L
        D
        Name
        Units
    end
    
    methods
        function C = Data2Cell(obj)
            C = DataConv(obj);
        end
    end
end

