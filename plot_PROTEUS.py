#!/bin/env python
# -*- coding: utf-8 -*-

#==============================================================================
# NAME : plot_PROTEUS
# AUTH : Guangchun Zhang
# EMAI : zhan2833@purdue.edu
#------------------------------------------------------------------------------
# VERS : 0.1
# TIME : 04/27/2017
# INFO : Plot lines of iteration information.
#==============================================================================

import os
import re
import sys
import optparse
import itertools
import matplotlib as mpl
# mpl.use('Qt5Agg')
import matplotlib.pyplot as plt

###############################################################################
def stripComment(s):
    return re.subn(r'!.*', '', s)[0]
###############################################################################

def parsePROTEUSInput(inp_file):
    '''
    Usage: parsePROTEUSInput(inp_file). inp_file is the name of PROTEUS input file.
    Function: Parse PROTEUS input file and return a dictionary in which KEY is the
              argument name and VALUE is the argumen value.
    '''
    input_arg_dic = {}
    with open(inp_file) as f:
        for line in f:
            line_list = stripComment(line).split()

            if len(line_list) > 2:
                key = line_list[0]
                value = line_list[1:]
            elif len(line_list) == 2:
                key   = line_list[0]
                value = line_list[1]
            elif len(line_list) == 1:
                sys.stderr.write("Missed argument value: %s"%key)
                sys.exit(1)
            else:
                continue

            input_arg_dic[key] = value

    return input_arg_dic

#==============================================================================
class IterationError():
    '''
    Class for holding iteraton message.
    '''
    #--------------------------------------------------------------------------
    def __init__(self):
        self.name = ''
        self.iter_dic = {}

    #--------------------------------------------------------------------------
    def setName(self, name):
        self.name = name

    #--------------------------------------------------------------------------
    def getName(self):
        return self.name

    #--------------------------------------------------------------------------
    def add(self, iter_key, iter_value):
        self.iter_dic[iter_key] = iter_value

    #--------------------------------------------------------------------------
    def getDict(self):
        return self.iter_dic

#==============================================================================
class IterationErrorList():
    '''
    class for holding list of IterationError.
    '''

    #--------------------------------------------------------------------------
    def __init__(self):
        self.iter_list = []

    #--------------------------------------------------------------------------
    def filter(self, solver_name, error_name):
        if len(self.iter_list) == 0:
            sys.stderr.write("Length of error list is zero. No value to filter.")
            sys.exit(1)

        # if not self.iter_list[0].getDict().has_key(error_name):
        if error_name not in self.iter_list[0].getDict():
            sys.stderr.write("%s can not be found!"%error_name)
            sys.exit(1)

        error_list = []

        for iter_error in self.iter_list:
            if (iter_error.getName() == solver_name):
                error_list.append(iter_error.getDict()[error_name])
                # print 'Name: ', iter_error.getName()

        return error_list
    #--------------------------------------------------------------------------
    def getError(self, error_name):
        if len(self.iter_list) == 0:
            sys.stderr.write("Length of error list is zero. No value to filter.\n")
            sys.exit(1)

        if not self.iter_list[0].getDict().has_key(error_name):
            sys.stderr.write("%s can not be found!\n"%error_name)
            sys.exit(1)

        error_list = []

        for iter_error in self.iter_list:
            error_list.append(iter_error.getDict()[error_name])

        return error_list

    #--------------------------------------------------------------------------
    def append(self, iter_error_dic):
        self.iter_list.append(iter_error_dic)

    #--------------------------------------------------------------------------
    def getIterNum(self):
        return len(self.iter_list)

    #--------------------------------------------------------------------------
    def getList(self):
        return self.iter_list

#==============================================================================
# Class for plot. Holding values for a single matplotlib.pyplot.plot()
# <Having not been used yet!>
class Error():
    def __init__(self, label='', error_name='', error_value=[]):
        self.label = label
        self.error_name = error_name
        self.error_value = error_value

    def getErrorName():
        return self.error_name

    def getErrorValue():
        return self.error_value

    def getLabel():
        return self.label

#==============================================================================
def parseMOCEXOutput(file_name):
    iter_error_list = IterationErrorList()

    BEGIN_SIGN = "BEGINNING OF EIGENVALUE SOLVE"
    END_SIGN   = "[ MOCEX ]"

    solver_name = ''

    f = open(file_name, 'r')
    offset = 0
    for line in f:
        if BEGIN_SIGN in line:
            break
        else:
            offset += len(line)
    f.seek(0)

    f.seek(offset)

    for i in range(3): f.readline()

    for line in f:
        error_dict = IterationError()
        line = line.rstrip()
        if END_SIGN in line: break
        line_list = line.split('|')
        # print line_list

        if 'MOCEX' in line_list[0] or 'LEGACY' in line_list[0]:
            solver_name = "MOC"
        elif 'CMFD' in line_list[0]:
            solver_name = 'CMFD'
        else:
            continue
            sys.stderr.write("Error occured when parsing this line:")
            sys.stderr.write(line)
            sys.exit(1)

        error_dict.setName(solver_name)
        # Total_Time_For_Iteration
        error_dict.add('Seconds', float(line_list[1]))
        # current iteration for MOCEX. total outer iteration numbers for CMFD
        error_dict.add('Itr', int(line_list[2]))
        # eigenvalue
        error_dict.add('Eigenvalue', float(line_list[3]))
        # Error_Eigenvalue
        error_dict.add('Error', float(line_list[4]))
        # Error_RMS_Fission
        error_dict.add('Fiss Err', float(line_list[5]))
        # Error_RMS_Flux
        error_dict.add('Flux Err', float(line_list[6]))
        # Estimated_DomRatio
        error_dict.add('Dom', float(line_list[7]))
        # Error_WGS_K_Maximum and Error_WGS_K_Group
        # Maximum of relative L2 norm error of flux and its group number.
        # Error_WGS_K_Maximum = MAX(ResidualNrom/RHSNorm)
        # Error_RMS_Flux = Error_WGS_K_Maximum
        error_dict.add('Error_WGS_K_Maximum', float(line_list[8].split()[0]))
        error_dict.add('Error_WGS_K_Group', int(line_list[8].split()[1]))
        # Iter_WGS_K_Cumulative. Total GMRES iteration count for an outer iteration.
        error_dict.add('Iter_WGS_K_Cumulative', int(line_list[9]))
        # Iter_WGS_K_Maximum and Iter_WGS_K_Maximum_Group
        error_dict.add('Iter_WGS_K_Maximum', int(line_list[10].split()[0]))
        error_dict.add('Iter_WGS_K_Maximum_Group', int(line_list[10].split()[1]))
        # Min WGS K iterations and its group
        error_dict.add('Iter_WGS_K_Minimum', int(line_list[11].split()[0]))
        error_dict.add('Iter_WGS_K_Minimum_Group', int(line_list[11].split()[1]))
        # Pseudo Error
        # if len(line_list) > 12:
        #     error_dict.add('PError', float(line_list[12].split()[0]))
        #     error_dict.add('PError_Group', int(line_list[12].split()[1]))

        iter_error_list.append(error_dict)

    return iter_error_list

#==============================================================================
def parseMOCEXResidual(file_name, outer=1, group=1):
    print('outer = ', outer, 'group = ', group)
    residual = []
    residual_list = []

    current_outer = 0
    current_group = 0

    BEGIN_SIGN = "BEGINNING OF EIGENVALUE SOLVE"
    END_SIGN   = "[ MOCEX ]"

    solver_name = ''

    f = open(file_name, 'r')
    offset = 0
    for line in f:
        if BEGIN_SIGN in line:
            break
        else:
            offset += len(line)
    f.seek(0)

    f.seek(offset)

    for i in range(3): f.readline()

    for line in f:
        line = line.rstrip()
        if END_SIGN in line: break
        # print line_list

        if 'CMFD' in line:
            current_outer = current_outer + 1
            continue

        if current_outer == outer:
            if 'Group' in line:
                linelist = line.split()
                current_group = current_group + 1
                assert(current_group == int(linelist[-1]))
                continue

        if current_outer == outer and current_group == group:
            if 'Right Hand Side Norm' in line:
                RHSNorm = float(line.split()[-1])
            elif 'residual' in line:
                residual_list.append(float(line.split()[-1])/RHSNorm)
            # iteration has converged
            elif 'convergence' in line:
                # residual_list.append(residual)
                break
            else:
                continue

    return residual_list
#==============================================================================
class MyParser(optparse.OptionParser):
    def format_epilog(self, formatter):
        return self.epilog

#==============================================================================
def parse_arguments():
    usage = 'Usage: %prog [options]'
    description = 'Plot iteration convergence behavior of MOCEX or other PROTEUS solvers.'
    epilog = '''
ErrorNameList:
  Seconds, Itr, Eigenvalue, Error, 'Fiss Err', 'Flux Err', Dom,
  Error_WGS_K_Maximum, Error_WGS_K_Group, Iter_WGS_K_Cumulative,
  Iter_WGS_K_Maximum, Iter_WGS_K_Maximum_Group, Iter_WGS_K_Minimum,
  Iter_WGS_K_Minimum_Group, PError, PError_Group;
  residual;

Warning:
  <-u> and <-g> must be used with <-r>.

Examples:
  %s --mocoutput mocex.out -s MOC -e Eigenvalue
  %s --mocoutput mocex.out -s MIX -e 'Fiss Err'
  %s --mocoutput "mocex_1.out mocex_2.out" -e Error
  %s --mocoutput "mocex_1.out mocex_2.out" -r -t 2 -g 1

'''%(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])
    parser = MyParser(usage=usage, description=description, epilog = epilog)

    parser.add_option('-o', '--output',
            action='store', dest='output',
            help='Output figure name. Default: outptu.', default='output')

    parser.add_option('-t', '--format',
            action='store', dest='format',
            help='Fomat of figure. Default: png.', default='png')

    parser.add_option('--mocinput',
            action='store', dest='mocinput',
            help='Input file of MOCEX. Default: INPUT. It is not used in this version.', default='INPUT')

    parser.add_option('--mocoutput',
            action='store', dest='mocoutput',
            help='Output file of MOCEX. May contain multiple output files. Default: OUTPUT.', default='OUTPUT')

    parser.add_option('-s', '--solver',
            action='store', dest='solver',
            help='Which solver error you want to plot [MOC|CMFD|MIX].\
                    Could be one or must be equal to mocoutput. Default: MOC', default='MOC')

    parser.add_option('-e', '--errorname',
            action='store', dest='errorname',
            help='Error name you want to plot. Default: Error.\
                    Could be one or must be equal to mocoutput. Separator: Semicolon (;)', default='Fiss Err')

    parser.add_option('-l', '--label',
            action='store', dest='label',
            help='Set label name of curve. Default: Output file name. Separator: Semicolon (;)')

    parser.add_option('-r', '--residual',
            action='store_true', dest='residual',
            help='Plot residual convergence of Group g and Outer n. Default: False', default=False)

    parser.add_option('-u', '--outer',
            action='store', dest='outer', type='int',
            help='Outer iteration index to plot the residual. Must be used with <-r>.', default=1)

    parser.add_option('-g', '--group',
            action='store', dest='group', type='int',
            help='Group index to plot the residual. Must be used with <-r>', default=1)

    #--------------------------------------------------------------------------
    # the following arguments should be used together.
    # parser.add_option('--outer',
    #         action='store', dest='outer',
    #         help='The Outer iteration number to be plotted. Should be used with -e=<Residual> option. Default:1', default=1)

    return parser.parse_args()

#==============================================================================
# def plot():
#==============================================================================

# Using lambda function to replace this ?
def getSolverName(solver_list, idx):
    if len(solver_list) == 1:
        return solver_list[0]
    else:
        return solver_list[idx]

def getErrorName(errorname_list, idx):
    if len(errorname_list) == 1:
        return errorname_list[0]
    else:
        return errorname_list[idx]

def getColorMap(iter_error_list_list, idx):
    color_map = []

    if iter_error_list_list[idx].getList()[0].getName() == 'MOC':
        color_map = ['cyan', 'red']
    elif iter_error_list_list[idx].getList()[0].getName() == 'CMFD':
        color_map = ['red', 'cyan']
    else:
        sys.stderr.write('Not implemted yet.\n')
        sys.exit(1)

    return color_map

#==============================================================================
def main():
    options, args = parse_arguments()

    # -------------------------------------------------------------------------
    # For now we don't use input file. So we don't check it.
    # mocinput_file_list = options.mocinput.split()
    # for mocinput_file in mocinput_file_list:
    #     if not os.path.isfile(options.mocinput):
    #         sys.stderr.write('MOCEX input file %s not existed!\n'%mocinput_file)
    #         sys.exit(1)

    # -------------------------------------------------------------------------
    mocoutput_file_list = options.mocoutput.split()
    for mocoutput_file in mocoutput_file_list:
        if not os.path.isfile(mocoutput_file):
            sys.stderr.write('MOCEX output file %s not existed!\n'%mocoutput_file)
            sys.exit(1)

    #--------------------------------------------------------------------------
    # get label for each curve. default: moc output file name
    label_list = []
    if options.label is None:
        label_list = mocoutput_file_list
    else:
        label_list = options.label.split(';')

    #--------------------------------------------------------------------------
    # solver list.
    # Either 1 or equal to the number of mocoutput.
    solver_list = options.solver.split()

    #--------------------------------------------------------------------------
    # error name list.
    # either 1 or equal to the number of mocoutput.
    errorname_list = options.errorname.split(';')

    # -------------------------------------------------------------------------

    iter_error_file = []
    error_list = []
    xdata = []
    xdata_len= []

    # ------------------------------------------------------------------------
    # parse each output file.
    if options.residual:
        for mocoutput_file in mocoutput_file_list:
            iter_error_file.append(parseMOCEXResidual(mocoutput_file, options.outer, options.group))
    else:
        for mocoutput_file in mocoutput_file_list:
            iter_error_file.append(parseMOCEXOutput(mocoutput_file))

    for idx in range(len(iter_error_file)):
        iter_error_list = iter_error_file[idx]

        solver_name = getSolverName(solver_list, idx)
        error_name  = getErrorName(errorname_list, idx)

        if not options.residual:
            if solver_name == 'MIX':
                error_list.append(iter_error_list.getError(error_name))
            else:
                error_list.append(iter_error_list.filter(solver_name, error_name))
        else:
            error_list = iter_error_file

    for idx in range(len(error_list)):
        print('file   : ', mocoutput_file_list[idx])
        print('solver : ', getSolverName(solver_list,idx))
        print('error  : ', options.errorname.split(';')[idx])
        print('length : ', len(error_list[idx]))
        print('value  :',  error_list[idx])
        print('sum    :',  sum(error_list[idx]))

    marker = itertools.cycle(( 'o', '^', '*', '<', '>', ',', '+', '.', 'v' ))

    # must be put before all string settings, like title, lable settings...
    font = {'family':'times', 'size':14}
    mpl.rc('font', **font)

    plt.figure(num=1, figsize=(9.708,6))
    ax = plt.gca()
    plt.subplot(111)
    plt.title('Converge info of %s'%options.errorname)
    plt.xlabel('Iteration Number')
    if not options.residual:
        plt.ylabel(options.errorname)
    else:
        plt.ylabel(r'$\Vert{r}\Vert_2/\Vert{b}\Vert_2$')

    error_name = options.errorname
    for idx in range(len(error_list)):
        error = error_list[idx]
        xdata = [i for i in range(1, len(error)+1)]
        xdata_len.append(len(xdata))
        solver_name = getSolverName(solver_list, idx)
        if error_name=='Eigenvalue' or error_name=='Dom' or error_name=='Iter_WGS_K_Cumulative' \
           or error_name=='Iter_WGS_K_Maximum' or error_name=='Iter_WGS_K_Maximum_Group' \
           or error_name=='Iter_WGS_K_Minimum' or error_name=='Iter_WGS_K_Minimum_Group':
            if solver_name != 'MIX':
                plt.plot(xdata, error, marker=next(marker), markersize=6,
                    label=label_list[idx])
            else:
                plt.scatter(xdata, error, s=80, color=getColorMap(iter_error_file, idx), marker=next(marker))
                plt.plot(xdata, error, label=label_list[idx])
        else:
            if solver_name != 'MIX':
                plt.semilogy(xdata, error, marker=next(marker), markersize=6,
                    label=label_list[idx])
            else:
                ax.set_yscale('log')
                plt.scatter(xdata, error, s=80, color=getColorMap(iter_error_file,idx), marker=next(marker))
                plt.semilogy(xdata, error, label=label_list[idx])

    plt.legend(loc='best')

    plt.xlim(0, max(xdata_len)+1)
    xticks = [i for i in range(max(xdata_len)+1)]
    # plt.ylim(0.9,1.15)
    # ax.set_xticks(xticks)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))

    plt.grid(b=True, which='major',  linestyle='--')
    plt.show()
    # plt.savefig(options.output+'.'+options.format, fmt=options.format)


#==============================================================================
if __name__ == '__main__':
    main()
