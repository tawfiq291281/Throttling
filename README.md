# Throttlinglatest
openfoamv2412
The throttling is a very interesting transient mass transfer phenomena in which due to high pressure difference a flowing liquid turns into vapor phase. This process is widely used in refrigeration and air conditioning system which operates on vapor compression refrigeration thermodynamic cycle. 

In Open-foam V-2412 there are different multiphase solvers which can capture phase change but in this study I opted to use a customized solver based on interphaseChangeFoam which is a VOF based multiphase solver designed for incompressible fluids. Though the solver belongs to incompressible family here it is used in conjunction with Peng-Robinson equation of state which helps to capture the compressibility of vapor phase.

But the main architecture of the solver depends on Lee’s phase change model. In which depending on the saturation temperature phase change of refrigerant occurs. And this saturation temperature depends on Antoine’s equation, which is a empirical corelation used to estimate the vapor pressure as a function of temperature. Here meshing is performed using Ansys fluent meshing platform and K-omega-SST turbulence model is used to capture turbulence. Results are displayed below where inlet temperature of R-32 is 313 K with mass flow rate of 0.025 kg/s.  
