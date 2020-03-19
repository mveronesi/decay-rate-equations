import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_acceptance import plot_acceptance
from declare_dec import dm, dg, gs, r, delta, gamma, beta
# acceptance
a_acc=widgets.BoundedFloatText(value=1.5, min=0, max=100, step=0.01, continuous_update=False, description=r'$a$')
n_acc=widgets.BoundedFloatText(value=1.5, min=0, max=100, step=0.01, continuous_update=False, description=r'$n$')
b_acc=widgets.BoundedFloatText(value=0.05, min=0, max=1, step=0.001, continuous_update=False, description=r'$b$')
beta_acc=widgets.BoundedFloatText(value=0.03, min=0, max=0.1, step=0.001, continuous_update=False, description=r'$\beta$')
cutoff_acc=widgets.BoundedFloatText(value=0.2, min=0, max=20, step=0.01, continuous_update=False, description=r'$t_{cut}$ [ps]')
# plotting
xmin=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc=widgets.BoundedFloatText(value=2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{dec}$ [a.u.]')
y_mix=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name=widgets.Text(value='plots/acceptance.eps',continuous_update=False)
save=widgets.ToggleButton(value=False,description='Save')
fold_amix = widgets.Checkbox(value=True,description=r'fold $A_{mix}$',disabled=False)
## Decay Time Acceptance
wbox_acc = HBox([VBox([a_acc,n_acc,b_acc,beta_acc,cutoff_acc]),VBox([xmin,xmax,y_osc,y_mix,fold_amix]),VBox([name,save])])
acc_pars = interactive_output(plot_acceptance,{'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                                'dm':dm,'dg':dg,'gs':gs,
                                                'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                                'xmin':xmin,'xmax':xmax,'y_osc':y_osc,'y_mix':y_mix,
                                                'name':name,'save':save,'fold_amix':fold_amix})
