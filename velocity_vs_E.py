import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def plot(x, y, y_err, xlabel: str, ylabel: str, title: str, filename):
    matplotlib.rcParams.update({'font.size': 24})
    plt.figure(figsize=(12, 8))
    plt.errorbar(x, y, y_err, fmt='o--', capsize=7, capthick=2, ecolor='red')
    plt.gca().invert_yaxis()
    plt.xscale('log')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.grid()
    # plt.savefig(f'output/{filename}', bbox_inches='tight')


def run_cmd(fields, name):
    # Хранилище результатов
    results = []
    
    for E in fields:
        print(f"Running E = {E} V/m...")
        
        # Запускаем программу
        cmd = f"wsl -d Ubuntu ./src/phosphorene 1000 300 {E[0]} {E[1]} 0 1e-11"
        output = subprocess.check_output(cmd, shell=True, text=True)
        
        # Парсим вывод
        data = {}
        data['Ex'] = E[0]
        data['Ey'] = E[1]
        
        for line in output.split('\n'):
            # Средняя скорость
            if "Average velocity:" in line:
                # Формат: "Average velocity: { vx, vy } μm/ps"
                parts = line.split('{')[1].split('}')[0].split(',')
                data['vx'] = float(parts[0].strip()) * 1e6
                data['vy'] = float(parts[1].strip()) * 1e6
            
            # Стандартное отклонение
            if "std:" in line:
                # Формат: "std: { stdx, stdy } μm/ps"
                parts = line.split('{')[1].split('}')[0].split(',')
                data['stdx'] = float(parts[0].strip()) * 1e6
                data['stdy'] = float(parts[1].strip()) * 1e6
                
            if "Scattering rates:" in line:
                # Формат: "Scattering rates: [sr0, ...] 1/s"
                parts = line.split('[')[1].split(']')[0].split(',')
                for part in range(0, len(parts)):
                    data[f'SR_{part}'] = float(parts[part].strip())
                    
        results.append(data)
    
    # Создаем DataFrame
    df = pd.DataFrame(results)
    
    # Сохраняем в CSV
    df.to_csv(f'output/simulation_results_{name}.csv', index=False)
    print("\nРезультаты сохранены в simulation_results.csv")
    print(df.head())
    
    return df

fields = [[1e2, 0], [1e3, 0], [1e4, 0]]
data = run_cmd(fields, 'Ex')

# data = pd.read_csv('output/simulation_results_Ex.csv')
plot(data.Ex, data.vx, data.stdx, '$E_x$, В/м', '$\\upsilon_x$, м/с', 'ВАХ', 'Vx(Ex).png')
plot(data.Ex, data.vx/data.Ex, data.stdx/data.Ex, '$E_x$, В/м', '$\\mu_x$, м$^2$/(В·с)', 'Подвижность', 'Mu(Ex).png')

# fields = [[0, 1e2], [0, 1e3], [0, 1e4], [0, 1e5], [0, 1e6], [0, 1e7]]
# data = run_cmd(fields, 'Ey')

# # data = pd.read_csv('output/simulation_results_Ey.csv')
# plot(data.Ey, data.vy, data.stdy, '$E_y$, В/м', '$\\upsilon_y$, м/с', 'ВАХ', 'Vy(Ey).png')
# plot(data.Ey, data.vy/data.Ey, data.stdy/data.Ey, '$E_y$, В/м', '$\\mu_y$, м$^2$/(В·с)', 'Подвижность', 'Mu(Ey).png')