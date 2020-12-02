import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_decrate_cpt import plot_decrate_cpt
# bmix
dm=widgets.BoundedFloatText(value=17.749, min=0, max=20, step=0.001, continuous_update=False, description=r'$\Delta m_{s}$ [ps$^{-1}$]')
dg=widgets.BoundedFloatText(value=0.083, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Delta \Gamma_{s}$ [ps$^{-1}$]')
gs=widgets.BoundedFloatText(value=0.6608, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Gamma_{s}$ [ps$^{-1}$]')
# afs- cpt
afs=widgets.BoundedFloatText(value=0, min=0, max=1, step=0.001, continuous_update=False, description=r'$a_{fs}$')
rez=widgets.BoundedFloatText(value=0, min=0, max=1, step=0.001, continuous_update=False, description=r'$Re(z)$')
imz=widgets.BoundedFloatText(value=0, min=0, max=1, step=0.001, continuous_update=False, description=r'$Im(z)$')
# plotting
xmin=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{dec}$ [a.u.]')
y_mix=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name=widgets.Text(value='plots/decrate_cpt.eps',continuous_update=False)#description='Save as:')
save=widgets.ToggleButton(value=False,description='Save')
# select
b_f = widgets.Checkbox(value=True,description=r'$B$$\to$$f$',disabled=False)
bbar_f = widgets.Checkbox(value=True,description=r'$\overline{B}$$\to$$f$',disabled=False)
bbar_fbar = widgets.Checkbox(value=True,description=r'$\overline{B}$$\to$$\overline{f}$',disabled=False)
b_fbar = widgets.Checkbox(value=True,description=r'$B$$\to$$\overline{f}$',disabled=False)
asymm =widgets.Dropdown(
    options=['MIX', 'CP', 'CPT'],
    value='MIX',
    description='Asymm',
    disabled=False,
)
## Decay Rate
wbox_cpt = HBox([VBox([gs,dg,dm]),VBox([rez,imz,afs,asymm]),VBox([xmin,xmax,y_osc,y_mix]),
                 VBox([b_f,bbar_f,bbar_fbar,b_fbar]),VBox([name,save])])
decrate_cpt = interactive_output(plot_decrate_cpt,{'dm':dm,'dg':dg,'gs':gs,
                                                'afs':afs,'rez':rez,'imz':imz,
                                                'xmin':xmin,'xmax':xmax,
                                                'y_osc':y_osc,'y_mix':y_mix,
                                                'name':name,'save':save,
                                                'b_f':b_f,
                                                'bbar_f':bbar_f,
                                                'bbar_fbar':bbar_fbar,
                                                'b_fbar':b_fbar,
                                                'asymm':asymm})
