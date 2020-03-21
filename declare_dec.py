import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_decrate import plot_decrate
# bmix - cp parameters
dm=widgets.BoundedFloatText(value=17.757, min=0, max=20, step=0.001, continuous_update=False, description=r'$\Delta m_{s}$ [ps$^{-1}$]')
dg=widgets.BoundedFloatText(value=0.085, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Delta \Gamma_{s}$ [ps$^{-1}$]')
gs=widgets.BoundedFloatText(value=0.664, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Gamma_{s}$ [ps$^{-1}$]')
r=widgets.BoundedFloatText(value=0, min=0, max=1, step=0.1, continuous_update=False, description=r'$r$')
delta=widgets.BoundedFloatText(value=0, min=0, max=360, step=1, continuous_update=False, description=r'$\delta$ [$\circ$]')
gamma=widgets.BoundedFloatText(value=0, min=0, max=360, step=1, continuous_update=False, description=r'$\gamma$ [$\circ$]')
beta=widgets.BoundedFloatText(value=0, min=-1000, max=1000, step=1, continuous_update=False, description=r'$\beta_{s}$ [mrad]')
# plotting
xmin=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc=widgets.BoundedFloatText(value=2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{dec}$ [a.u.]')
y_mix=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name=widgets.Text(value='plots/decrate.eps',continuous_update=False)#description='Save as:')
save=widgets.ToggleButton(value=False,description='Save')
# select
b_f = widgets.Checkbox(value=True,description=r'$B$$\to$$f$',disabled=False)
bbar_f = widgets.Checkbox(value=True,description=r'$\overline{B}$$\to$$f$',disabled=False)
bbar_fbar = widgets.Checkbox(value=True,description=r'$\overline{B}$$\to$$\overline{f}$',disabled=False)
b_fbar = widgets.Checkbox(value=True,description=r'$B$$\to$$\overline{f}$',disabled=False)
fold_amix = widgets.Checkbox(value=True,description=r'fold $A_{mix}$',disabled=False)
## Decay Rate
wbox_dec = HBox([VBox([gs,dg,dm,fold_amix]),VBox([r,delta,gamma,beta]),VBox([xmin,xmax,y_osc,y_mix]),
                 VBox([b_f,bbar_f,bbar_fbar,b_fbar]),VBox([name,save])])
decrate_pars = interactive_output(plot_decrate,{'dm':dm,'dg':dg,'gs':gs,
                                                'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                                'xmin':xmin,'xmax':xmax,'y_osc':y_osc,'y_mix':y_mix,
                                                'name':name,'save':save,
                                                'b_f':b_f,
                                                'bbar_f':bbar_f,
                                                'bbar_fbar':bbar_fbar,
                                                'b_fbar':b_fbar,
                                                'fold_amix':fold_amix})
