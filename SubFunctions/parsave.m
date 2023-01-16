function parsave(fname,data)

%var_name=genvarname(inputname(2)); 
var_name = matlab.lang.makeValidName(inputname(2)) % inputname 
inputname(2)
eval([var_name '=data'])

try 
    save(fname,var_name,'-append') 
catch 
    save(fname,var_name) 
end

