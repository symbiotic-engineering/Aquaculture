import pandas as pd
import numpy as np

def load_wave_data(file_name):
    if not file_name:
        wave_period = 8.33 * np.ones(8760)
        wave_height = 1.40 * np.ones(8760)
        return wave_period , wave_height
    else:
        df = pd.read_csv(file_name)
        wave_period = df['Peak Period']
        wave_height = df['Significant Wave Height']
        return np.array(wave_period.values), np.array(wave_height.values)
    
    
