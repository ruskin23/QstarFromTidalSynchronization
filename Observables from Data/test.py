

import scipy

mass_array = scipy.linspace(1, 20, 20)





t_array = scipy.empty(len(mass_array))

for m_index, m_value in enumerate(mass_array):
    t_array[m_index] = m_index*10

"""t_array[:] = [x - 10 for x in t_array]"""


print(t_array)





