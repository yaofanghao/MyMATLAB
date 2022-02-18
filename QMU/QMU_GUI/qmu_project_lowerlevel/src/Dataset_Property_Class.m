classdef Dataset_Property_Class
    properties 
        Distribution
        Distribution_value
    end
    methods
        function obj = Dataset_Property_Class(n)
            if ( n == 1)
                obj.Distribution = 'Gaussian distribution or Logarithmic Gaussian distribution';
                obj.Distribution_value = 1;
            else
                obj.Distribution = 'Logistic distribution or Logarithmic Logistic distribution';
                obj.Distribution_value = 0;
            end
        end
    end
end
            
        