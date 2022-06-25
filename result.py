import numpy as np 
import matplotlib.pyplot as plt
from scipy import optimize as opt

gen = 150
x = np.arange(0, gen)
y = []

#basic settings
plt.xlabel('Generarion')
plt.ylabel('Max fitness')

#plt.plot(x, y, 'r-', label="Roulette Select") # - line o round 
plt.plot(x, y, 'r-', label="RouletteSelect") # - line o round  
"""
plt.plot(x, y0, 'r-', label="RouletteSelect") # - line o round   
plt.plot(x, y1, 'g-', label="TournamentSelect") # - line o round   
plt.plot(x, y2, 'b-', label="ExpectedValueSelect") # - line o round   
plt.plot(x, y3, 'c-', label="RankingSelect (Proportional)") # - line o round   
plt.plot(x, y4, 'm-', label="RankingSelect (Inverse proportional)") # - line o round   
plt.plot(x, y5, 'y-', label="RankingSelect (Quadratic function)") # - line o round 
"""
"""
plt.plot(x, y6, 'k-', label="num=7") # - line o round   
plt.plot(x, y7, color='#377eb8', label="num=8") # - line o round   
plt.plot(x, y8, color='#e41a1c', label="num=9") # - line o round   
plt.plot(x, y9, color='#ff7f00', label="num=10") # - line o round  
""" 
#plt.scatter(gen_x, gen_y) # - line o round 
plt.legend(loc="best", fontsize=9) #if use label (loc = upper right)  
plt.grid(True)
plt.show()


