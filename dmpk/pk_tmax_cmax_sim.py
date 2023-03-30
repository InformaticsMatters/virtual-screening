# Copyright 2023 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Estimation of absorption from t1/2 and Tmax after po administration (one-compartment)
# Based on original work by Amit Kumar Garg <a.garg@sygnaturediscovery.com>


import argparse, collections
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import utils


def generatePlot(t_hf, t_hf_a, D, AUC, tn, quiet=False, plot_height=4, plot_width=10, font_size=12, basename='output'):

    utils.expand_path(basename + '.png')

    kel = math.log(2)/t_hf
    ka = math.log(2)/ t_hf_a
    Tmax = (math.log(ka)-math.log(kel))/(ka-kel)
    Cmax = math.exp(-kel*Tmax)*kel*AUC
    V_F = D/kel/AUC

    outputs = collections.OrderedDict()
    outputs['Tmax(hr)'] = utils.round_to_significant_number(Tmax, 3)
    outputs['Cmax(mg/L)'] = utils.round_to_significant_number(Cmax, 3)
    outputs['Kel(hr-1)'] = utils.round_to_significant_number(kel, 3)
    outputs['Ka(hr-1)'] = utils.round_to_significant_number(ka, 3)
    outputs['V/F(L)'] = utils.round_to_significant_number(V_F, 3)

    if not quiet:
        utils.log('------------------------------------------------------------------------------------------')
        utils.log('kel \t', kel)
        utils.log('ka \t', ka)
        utils.log('Tmax \t', Tmax)
        utils.log('Cmax \t', Cmax)
        utils.log('V_F \t', V_F)
        utils.log('------------------------------------------------------------------------------------------')

    b_time = []
    c_cp = []
    d_perc = []
    for i in range(0, 101):
        if i == 0:
            b_time.append(0)
        else:
            b_time.append(b_time[i-1] + tn/100)

        c_cp.append(ka*D/V_F/(ka-kel)*(math.exp(-kel*b_time[i])-math.exp(-ka*b_time[i])))

        d_perc.append(100-100*math.exp(-ka*b_time[i]))

    # Creating the visualisation
    plt.figure(figsize=(plot_width, plot_height))
    plt.subplot(1, 2, 1)

    plt.plot(b_time, c_cp, linewidth=2, linestyle='dashed', color='coral')  # Plotting the observed data
    plt.xlabel('Time (h)', fontsize=font_size)
    plt.ylabel('Cp(mg/L)', fontsize=font_size)
    plt.title('cp Vs Time', color='coral', fontsize=font_size)
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(b_time,d_perc, linewidth=2, linestyle='dashed')  # Plotting the observed data
    plt.xlabel('Time (h)', fontsize=font_size)
    plt.ylabel('% Absorbed', fontsize=font_size)
    plt.title('%Absorbed Vs Time', color='dodgerblue', fontsize=font_size)
    plt.grid(True)

    # Fine-tune figure; make subplots farther from each other.
    # refine layout to better support different sizes

    plt.savefig(basename + '.png')
    with open(basename + '.txt', 'wt') as writer:
        writer.write('Inputs:\n')
        writer.write('  half-life = {}\n  absorption = {}\n  dose = {}\n  auc = {}\n  time = {}\n'.
                     format(t_hf, t_hf_a, D, AUC, tn))
        writer.write('Outputs:\n')
        for k, v in outputs.items():
            writer.write('  {} = {}\n'.format(k, v))

    return outputs


def main():

    # example usage:
    #   python -m dmpk.pk_tmax_cmax_sim --half-life 0.79 --absorption 0.5 --dose 0.14 --auc 0.88 --time 8

    parser = argparse.ArgumentParser(description='Tmax/Cmax simulation')
    parser.add_argument('--half-life', type=float, required=True, help='half life (hours)')
    parser.add_argument('--absorption', type=float, required=True, help='half life absorption (hours)')
    parser.add_argument('--dose', type=float, required=True, help='initial dose (mg)')
    parser.add_argument('--auc', type=float, required=True, help='AUC (mg/L*hr)')
    parser.add_argument('--time', type=float, required=True, help='time (h)')

    parser.add_argument('--plot-height', type=int, default=4, help='plot height (in)')
    parser.add_argument('--plot-width', type=int, default=10, help='plot width (in)')
    parser.add_argument('--font-size', type=int, default=12, help='font size (points)')

    parser.add_argument('-o', '--output', type=str, default='output', help='output file base name')

    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')

    args = parser.parse_args()
    utils.log("Tmax/Cmax simulation Args: ", args)

    ### execute #######################################################

    outputs = generatePlot(args.half_life, args.absorption, args.dose, args.auc, args.time,
                           plot_width=args.plot_width, plot_height=args.plot_height, font_size=args.font_size,
                           basename=args.output)


if __name__ == "__main__":
    main()
