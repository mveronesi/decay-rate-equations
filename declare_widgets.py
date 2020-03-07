import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_decrate import plot_decrate
from plot_acceptance import plot_acceptance
from plot_resolution import plot_resolution
from plot_tagging import plot_tagging
from plot_asymmetries import plot_asymmetries

# Widgets
# decay rate
dm=widgets.BoundedFloatText(value=17.757, min=0, max=20, step=0.001, continuous_update=False, description=r'$\Delta m_{s}$ [ps$^{-1}$]')
dg=widgets.BoundedFloatText(value=0.085, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Delta \Gamma_{s}$ [ps$^{-1}$]')
gs=widgets.BoundedFloatText(value=0.664, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Gamma_{s}$ [ps$^{-1}$]')
r=widgets.BoundedFloatText(value=0, min=0, max=1, step=0.1, continuous_update=False, description=r'$r$')
delta=widgets.BoundedFloatText(value=0, min=0, max=360, step=1, continuous_update=False, description=r'$\delta$ [$\circ$]')
gamma=widgets.BoundedFloatText(value=0, min=0, max=180, step=1, continuous_update=False, description=r'$\gamma$ [$\circ$]')
beta=widgets.BoundedFloatText(value=0, min=-100, max=100, step=1, continuous_update=False, description=r'$\beta_{s}$ [mrad]')
# acceptance
a_acc=widgets.BoundedFloatText(value=1.5, min=0, max=100, step=0.01, continuous_update=False, description=r'$a$')
n_acc=widgets.BoundedFloatText(value=1.5, min=0, max=100, step=0.01, continuous_update=False, description=r'$n$')
b_acc=widgets.BoundedFloatText(value=0.05, min=0, max=1, step=0.001, continuous_update=False, description=r'$b$')
beta_acc=widgets.BoundedFloatText(value=0.03, min=0, max=0.1, step=0.001, continuous_update=False, description=r'$\beta$')
cutoff_acc=widgets.BoundedFloatText(value=0.2, min=0, max=20, step=0.01, continuous_update=False, description=r'$t_{cut}$ [ps]')
# resolution
sigma_t=widgets.BoundedFloatText(value=0.045, min=0, max=1, step=0.001, continuous_update=False, description=r'$\sigma_{t}$ [ps]')
# tagging
omega=widgets.BoundedFloatText(value=0.35, min=0, max=0.5, step=0.01, continuous_update=False, description=r'$\omega_{tag}$')
d_omega=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\omega_{tag}$')
eff=widgets.BoundedFloatText(value=0.8, min=0, max=1, step=0.01, continuous_update=False, description=r'$\varepsilon_{tag}$')
d_eff=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\varepsilon_{tag}$')
# asymmetries
a_prod=widgets.BoundedFloatText(value=0, min=-1, max=1, step=0.01, continuous_update=False, description=r'$a_{prod}$')
a_det=widgets.BoundedFloatText(value=0, min=-1, max=1, step=0.01, continuous_update=False, description=r'$a_{det}$')
# plotting
## dec_rate
xmin_dec=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps$^{-1}$]')
xmax_dec=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps$^{-1}$]')
y_osc_dec=widgets.BoundedFloatText(value=2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{osc}$ [a.u.]')
y_mix_dec=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_dec=widgets.Text(value='plots/decrate.eps',continuous_update=False,description='Save as:')
save_dec=widgets.ToggleButton(value=False,description='Save')
## acceptance
xmin_acc=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps$^{-1}$]')
xmax_acc=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps$^{-1}$]')
y_osc_acc=widgets.BoundedFloatText(value=2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{osc}$ [a.u.]')
y_mix_acc=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_acc=widgets.Text(value='plots/acceptance.eps',continuous_update=False,description='Save as:')
save_acc=widgets.ToggleButton(value=False,description='Save')
## resolution
xmin_res=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps$^{-1}$]')
xmax_res=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps$^{-1}$]')
y_osc_res=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{osc}$ [a.u.]')
y_mix_res=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_res=widgets.Text(value='plots/resolution.eps',continuous_update=False,description='Save as:')
save_res=widgets.ToggleButton(value=False,description='Save')
## tagging
xmin_tag=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps$^{-1}$]')
xmax_tag=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps$^{-1}$]')
y_tag=widgets.BoundedFloatText(value=0.4, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{tag}$ [a.u.]')
y_untag=widgets.BoundedFloatText(value=0.2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{untag}$ [a.u.]')
y_mix_tag=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_tag=widgets.Text(value='plots/tagging.eps',continuous_update=False,description='Save as:')
save_tag=widgets.ToggleButton(value=False,description='Save')
## asymmetries
xmin_asymm=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps$^{-1}$]')
xmax_asymm=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps$^{-1}$]')
y_tag_asymm=widgets.BoundedFloatText(value=0.4, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{tag}$ [a.u.]')
y_untag_asymm=widgets.BoundedFloatText(value=0.2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{untag}$ [a.u.]')
y_mix_asymm=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_asymm=widgets.Text(value='plots/asymmetries.eps',continuous_update=False,description='Save as:')
save_asymm=widgets.ToggleButton(value=False,description='Save')

# Interactive boxes
## Decay Rate
wbox_dec = HBox([VBox([dm,dg,gs]),VBox([r,delta,gamma,beta]),VBox([xmin_dec,xmax_dec,y_osc_dec,y_mix_dec]),VBox([name_dec]),VBox([save_dec])])
decrate_pars = interactive_output(plot_decrate,{'dm':dm,'dg':dg,'gs':gs,
                                                'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                                'xmin':xmin_dec,'xmax':xmax_dec,'y_osc':y_osc_dec,'y_mix':y_mix_dec,
                                                'name':name_dec,'save':save_dec})
## Decay Time Acceptance
wbox_acc = HBox([VBox([a_acc,n_acc,b_acc,beta_acc,cutoff_acc]),VBox([xmin_acc,xmax_acc,y_osc_acc,y_mix_acc]),VBox([name_acc]),VBox([save_acc])])
acc_pars = interactive_output(plot_acceptance,{'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                                'dm':dm,'dg':dg,'gs':gs,
                                                'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                                'xmin':xmin_acc,'xmax':xmax_acc,'y_osc':y_osc_acc,'y_mix':y_mix_acc,
                                                'name':name_acc,'save':save_acc})
## Decay Time Resolution
wbox_res = HBox([VBox([sigma_t]),VBox([xmin_res,xmax_res,y_osc_res,y_mix_res]),VBox([name_res]),VBox([save_res])])
res_pars = interactive_output(plot_resolution,{'sigma_t':sigma_t,
                                                'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                                'dm':dm,'dg':dg,'gs':gs,
                                                'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                                'xmin':xmin_res,'xmax':xmax_res,'y_osc':y_osc_res,'y_mix':y_mix_res,
                                                'name':name_res,'save':save_res})
## Flavoure Tagging
wbox_tag = HBox([VBox([omega,d_omega,eff,d_eff]),VBox([xmin_tag,xmax_tag,y_tag,y_untag,y_mix_tag]),VBox([name_tag]),VBox([save_tag])])
tag_pars = interactive_output(plot_tagging,{'omega':omega,'d_omega':d_omega,'eff':eff,'d_eff':d_eff,
                                            'sigma_t':sigma_t,
                                            'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                            'dm':dm,'dg':dg,'gs':gs,
                                            'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                            'xmin':xmin_tag,'xmax':xmax_tag,'y_tag':y_tag,'y_untag':y_untag,'y_mix':y_mix_tag,
                                            'name':name_tag,'save':save_tag})
## Asymmetries
wbox_asymm = HBox([VBox([a_prod,a_det]),VBox([xmin_asymm,xmax_asymm,y_tag_asymm,y_untag_asymm,y_mix_asymm]),VBox([name_asymm]),VBox([save_asymm])])
asymm_pars = interactive_output(plot_asymmetries,{'a_prod':a_prod,'a_det':a_det,
                                            'omega':omega,'d_omega':d_omega,'eff':eff,'d_eff':d_eff,
                                            'sigma_t':sigma_t,
                                            'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                            'dm':dm,'dg':dg,'gs':gs,
                                            'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                            'xmin':xmin_asymm,'xmax':xmax_asymm,'y_tag':y_tag_asymm,'y_untag':y_untag_asymm,'y_mix':y_mix_asymm,
                                            'name':name_asymm,'save':save_asymm})