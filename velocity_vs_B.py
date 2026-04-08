import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def run_cmd(fields, E, name):
    # Хранилище результатов
    results = []
    
    for B in fields:
        print(f"Running B = {B} T...")
        
        # Запускаем программу
        cmd = f"wsl -d Ubuntu ./src/phosphorene 1000 300 {E} {B} 1e-11"
        output = subprocess.check_output(cmd, shell=True, text=True)
        
        # Парсим вывод
        data = {}
        data['B'] = B
        
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
                    
        results.append(data)
    
    # Создаем DataFrame
    df = pd.DataFrame(results)
    
    # Сохраняем в CSV
    df.to_csv(f'output/simulation_results_{name}.csv', index=False)
    print("\nРезультаты сохранены в simulation_results.csv")
    print(df.head())
    
    return df

def plot(x, y, y_err, xlabel: str, ylabel: str, title: str):
    matplotlib.rcParams.update({'font.size': 24})
    plt.figure(figsize=(12, 8))
    plt.errorbar(x, y, y_err, fmt='o--', capsize=7, capthick=2, ecolor='red')
    # plt.gca().invert_yaxis()
    plt.xscale('log')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.grid()
    plt.savefig('output/V(B).png', bbox_inches='tight')

fields = [1e-3, 1e-2, 1e-1, 1, 1e1]
E = [1e3, 0]
data = run_cmd(fields, E, 'B')

# data = pd.read_csv('output/simulation_results_B.csv')
plot(data.B, data.vy, data.stdy, '$B$, Т', '$\\upsilon_y$, м/с', 'эффект Холла')

