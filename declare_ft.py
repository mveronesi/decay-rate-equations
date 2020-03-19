import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_resolution import plot_resolution
from declare_dec import dm, dg, gs, r, delta, gamma, beta
from declare_acc import a_acc, n_acc, b_acc, beta_acc, cutoff_acc
from declare_res import sigma_t

# tagging
omega=widgets.BoundedFloatText(value=0.35, min=0, max=0.5, step=0.01, continuous_update=False, description=r'$\omega_{tag}$')
d_omega=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\omega_{tag}$')
eff=widgets.BoundedFloatText(value=0.8, min=0, max=1, step=0.01, continuous_update=False, description=r'$\varepsilon_{tag}$')
d_eff=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\varepsilon_{tag}$')
# plotting
xmin=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{dec}$ [a.u.]')
y_mix=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name=widgets.Text(value='plots/tagging.eps',continuous_update=False)
save=widgets.ToggleButton(value=False,description='Save')
# select
qt_plus = widgets.Checkbox(value=True,description=r'$q_t$=+1',disabled=False)
qt_min = widgets.Checkbox(value=True,description=r'$q_t$=-1',disabled=False)
qf_plus = widgets.Checkbox(value=True,description=r'$q_f$=+1',disabled=False)
qf_min = widgets.Checkbox(value=True,description=r'$q_f$=-1',disabled=False)
fold_amix = widgets.Checkbox(value=True,description=r'fold $A_{mix}$',disabled=False)
k_acc= widgets.Checkbox(value=True,description=r'Acceptance',disabled=False)
k_res= widgets.Checkbox(value=True,description=r'Resolution',disabled=False)

## Flavoure Tagging
wbox_tag = HBox([VBox([omega,d_omega,eff,d_eff,k_acc_tag,k_res_tag]),VBox([xmin_tag,xmax_tag,y_tag,y_untag,y_mix_tag]),VBox([name_tag]),VBox([save_tag])])
tag_pars = interactive_output(plot_tagging,{'omega':omega,'d_omega':d_omega,'eff':eff,'d_eff':d_eff,
                                            'sigma_t':sigma_t,
                                            'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                            'dm':dm,'dg':dg,'gs':gs,
                                            'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                            'xmin':xmin_tag,'xmax':xmax_tag,'y_tag':y_tag,'y_untag':y_untag,'y_mix':y_mix_tag,
                                            'k_acc':k_acc_tag,'k_res':k_res_tag,
                                            'name':name_tag,'save':save_tag})
