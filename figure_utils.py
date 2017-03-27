
def get_pretty_name(line):
    if line[0]=='m':
        return 'Ara-%s' % line[1]
    else:
        return 'Ara+%s' % line[1]

line_color_map = {'m5': '#4A1486', 'p2': '#807DBA', 'p4': '#084594','p1': '#4292C6', 'p5': '#005A32', 'm6': '#41AB5D', 'm1': '#8C2D04', 'm2': '#CC4C02', 'm3': '#B10026', 'm4': '#E31A16', 'p3': '#FC4E2A', 'p6': '#FD8D3C'}

def get_line_color(population):
    return line_color_map[population]

nonmutator_group_color = get_line_color('p4')
mutator_group_color = get_line_color('m3')

nonmutator_group_label = 'Nonmutators'
mutator_group_label = 'Mutators'

time_axis_label = "Generations"

time_xticks = [10000*i for i in xrange(0,7)]
time_xticklabels = ['%dk' % (10*i) for i in xrange(0,7)]
time_xticklabels[0] = '0'


var_type_color_map = {'synonymous': '#e6ab02',
                      'missense': '#66a61e',
                      'nonsense': '#d95f02',
                      'noncoding': '#7570b3',
                      'indel': '#a6761d',
                      'sv': '#1b9e77'}
                      

def get_var_type_color(var_type):
    return var_type_color_map[var_type]



def get_panel_label(c):
    
    return c.lower()
    #return c.upper()