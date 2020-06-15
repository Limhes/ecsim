#!/bin/sh

# Er:
#./equation_gen.sh "i_p = 0.4463 n F A C \sqrt{\frac{n F \nu D_{red}}{R T}}" formula_Er_ip.png
#./equation_gen.sh "E_p = \left( 1.109 + \ln{\sqrt{\frac{D_{red}}{D_{ox}}}} \right) \frac{R T}{n F}" formula_Er_Ep.png

# Ei:
#./equation_gen.sh "i_p = 0.4958 F A C \sqrt{\frac{(1 - \alpha) F \nu D_{red}}{R T}}" formula_Ei_ip.png
#./equation_gen.sh "E_p = \left( 0.780 - \ln{k_e} + \ln{\sqrt{ \frac{D_{red} \nu (1 - \alpha) F}{R T} }} \right) \frac{R T}{(1 - \alpha) F}" formula_Ei_Ep.png

# ErCr:
#./equation_gen.sh "i_p = 0.4463 F A C \sqrt{\frac{F \nu D}{R T}}" formula_ErCr_ip.png
#./equation_gen.sh "E_p = \left( 1.109 - \ln{\left( 1 + \frac{k_f}{k_b} \right)} \right) \frac{R T}{F}" formula_ErCr_Ep.png

# ErCi:
#./equation_gen.sh "i_p = 0.4958 F A C \sqrt{\frac{F \nu D}{R T}}" formula_ErCi_ip.png
#./equation_gen.sh "E_p = \left( 0.780 - 0.5 \ln{k_f} - 0.5 \ln{\sqrt{\frac{\nu F}{R T}} } \right) \frac{R T}{F}" formula_ErCi_Ep.png

# CE in DE:
#./equation_gen.sh "i_p = 0.4463 n F A C \sqrt{\frac{n F \nu D_{red}}{R T}}" formula_CE_DE_ip.png
#./equation_gen.sh "E_p = \left( 1.109 - \ln{ \frac{k_f}{k_f + k_b} } \right) \frac{R T}{F}" formula_CE_DE_Ep.png

# CE in KP:
./equation_gen.sh "i_{plateau} = F A C \sqrt{D k_b} \cdot k_f/k_b" formula_CE_KP_ip.png


rm -f formula.aux formula.log formula.pdf texput.log
