import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_resolution import plot_resolution
from declare_dec import dm, dg, gs, r, delta, gamma, beta
from declare_acc import a_acc, n_acc, b_acc, beta_acc, cutoff_acc

# resolution
sigma_t=widgets.BoundedFloatText(value=0.045, min=0, max=1, step=0.001, continuous_update=False, description=r'$\sigma_{t}$ [ps]')
# plotting
xmin=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{dec}$ [a.u.]')
y_mix=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name=widgets.Text(value='plots/resolution.eps',continuous_update=False)
save=widgets.ToggleButton(value=False,description='Save')
# select
qt_plus = widgets.Checkbox(value=True,description=r'$q_t$=+1',disabled=False)
qt_min = widgets.Checkbox(value=True,description=r'$q_t$=-1',disabled=False)
qf_plus = widgets.Checkbox(value=True,description=r'$q_f$=+1',disabled=False)
qf_min = widgets.Checkbox(value=True,description=r'$q_f$=-1',disabled=False)
fold_amix = widgets.Checkbox(value=True,description=r'fold $A_{mix}$',disabled=False)
k_acc= widgets.Checkbox(value=True,description=r'Acceptance',disabled=False)

## Decay Time Resolution
wbox_res = HBox([VBox([sigma_t,k_acc,fold_amix]),VBox([xmin,xmax,y_osc,y_mix]),
                 VBox([qt_plus,qt_min,qf_plus,qf_min]),
                 VBox([name]),VBox([save])])
res_pars = interactive_output(plot_resolution,{'sigma_t':sigma_t,
                                                'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                                'dm':dm,'dg':dg,'gs':gs,
                                                'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                                'xmin':xmin,'xmax':xmax,'y_osc':y_osc,'y_mix':y_mix,
                                                'k_acc':k_acc,
                                                'name':name,'save':save,
                                                'qt_plus':qt_plus,
                                                'qt_min':qt_min,
                                                'qf_plus':qf_plus,
                                                'qf_min':qf_min,
                                                'fold_amix':fold_amix})
