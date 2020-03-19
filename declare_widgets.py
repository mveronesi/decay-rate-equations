import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_decrate import plot_decrate
from plot_acceptance import plot_acceptance
from plot_resolution import plot_resolution
from plot_tagging import plot_tagging
from plot_asymmetries import plot_asymmetries
from plot_coeffs import plot_coeffs

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
k_acc_res= widgets.Checkbox(value=True,description=r'Acceptance',disabled=False)
# tagging
omega=widgets.BoundedFloatText(value=0.35, min=0, max=0.5, step=0.01, continuous_update=False, description=r'$\omega_{tag}$')
d_omega=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\omega_{tag}$')
eff=widgets.BoundedFloatText(value=0.8, min=0, max=1, step=0.01, continuous_update=False, description=r'$\varepsilon_{tag}$')
d_eff=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\varepsilon_{tag}$')
k_acc_tag= widgets.Checkbox(value=True,description=r'Acceptance',disabled=False)
k_res_tag= widgets.Checkbox(value=True,description=r'Resolution',disabled=False)
# asymmetries
a_prod=widgets.BoundedFloatText(value=-0.01, min=-1, max=1, step=0.001, continuous_update=False, description=r'$a_{prod}$')
a_det=widgets.BoundedFloatText(value=0.01, min=-1, max=1, step=0.001, continuous_update=False, description=r'$a_{det}$')
# decay rate coefficients
# qt = widgets.RadioButtons(options=[+1,0,-1],value=+1,description=r'$q_t$',disabled=False)
# qf = widgets.RadioButtons(options=[+1,-1],value=+1,description=r'$q_f$',disabled=False)
qt_plus = widgets.Checkbox(value=True,description=r'$q_t$=+1',disabled=False)
qt_min = widgets.Checkbox(value=False,description=r'$q_t$=-1',disabled=False)
qt_nul = widgets.Checkbox(value=False,description=r'$q_t$=0',disabled=False)
qf_plus = widgets.Checkbox(value=True,description=r'$q_f$=+1',disabled=False)
qf_min = widgets.Checkbox(value=False,description=r'$q_f$=-1',disabled=False)
k_dec = widgets.Checkbox(value=True,description=r'Decay',disabled=False)
k_acc= widgets.Checkbox(value=False,description=r'Acceptance',disabled=False)
k_res= widgets.Checkbox(value=False,description=r'Resolution',disabled=False)
k_tag= widgets.Checkbox(value=False,description=r'Tagging',disabled=False)
# k_prod= widgets.Checkbox(value=False,description=r'Production',disabled=False)
# k_det= widgets.Checkbox(value=False,description=r'Detection',disabled=False)
k_asymm= widgets.Checkbox(value=False,description=r'Asymmetries',disabled=False)
# plotting
## dec_rate
xmin_dec=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax_dec=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc_dec=widgets.BoundedFloatText(value=2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{osc}$ [a.u.]')
y_mix_dec=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_dec=widgets.Text(value='plots/decrate.eps',continuous_update=False,description='Save as:')
save_dec=widgets.ToggleButton(value=False,description='Save')
## acceptance
xmin_acc=widgets.BoundedFloatText(value=0, min=0, max=100, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax_acc=widgets.BoundedFloatText(value=5, min=0, max=100, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc_acc=widgets.BoundedFloatText(value=2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{osc}$ [a.u.]')
y_mix_acc=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_acc=widgets.Text(value='plots/acceptance.eps',continuous_update=False,description='Save as:')
save_acc=widgets.ToggleButton(value=False,description='Save')
## resolution
xmin_res=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax_res=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_osc_res=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{osc}$ [a.u.]')
y_mix_res=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_res=widgets.Text(value='plots/resolution.eps',continuous_update=False,description='Save as:')
save_res=widgets.ToggleButton(value=False,description='Save')
## tagging
xmin_tag=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax_tag=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_tag=widgets.BoundedFloatText(value=0.4, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{tag}$ [a.u.]')
y_untag=widgets.BoundedFloatText(value=0.2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{untag}$ [a.u.]')
y_mix_tag=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_tag=widgets.Text(value='plots/tagging.eps',continuous_update=False,description='Save as:')
save_tag=widgets.ToggleButton(value=False,description='Save')
## asymmetries
xmin_asymm=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax_asymm=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_tag_asymm=widgets.BoundedFloatText(value=0.4, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{tag}$ [a.u.]')
y_untag_asymm=widgets.BoundedFloatText(value=0.2, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{untag}$ [a.u.]')
y_mix_asymm=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{mix}$ [a.u.]')
name_asymm=widgets.Text(value='plots/asymmetries.eps',continuous_update=False,description='Save as:')
save_asymm=widgets.ToggleButton(value=False,description='Save')
## decay rate coefficients
xmin_coeffs=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax_coeffs=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_cosh_coeffs=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\cosh}$ [a.u.]')
y_sinh_coeffs=widgets.BoundedFloatText(value=0.02, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\sinh}$ [a.u.]')
y_cos_coeffs=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\cos}$ [a.u.]')
y_sin_coeffs=widgets.BoundedFloatText(value=1, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\sin}$ [a.u.]')
name_coeffs=widgets.Text(value='plots/coeffs.eps',continuous_update=False,description='Save as:')
save_coeffs=widgets.ToggleButton(value=False,description='Save')

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
wbox_res = HBox([VBox([sigma_t,k_acc_res]),VBox([xmin_res,xmax_res,y_osc_res,y_mix_res]),VBox([name_res]),VBox([save_res])])
res_pars = interactive_output(plot_resolution,{'sigma_t':sigma_t,
                                                'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                                                'dm':dm,'dg':dg,'gs':gs,
                                                'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                                                'xmin':xmin_res,'xmax':xmax_res,'y_osc':y_osc_res,'y_mix':y_mix_res,
                                                'k_acc':k_acc_res,
                                                'name':name_res,'save':save_res})
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

## Decay Rate Coefficients
wbox_coeffs_test = HBox([VBox([qt_plus,qt_nul,qt_min,qf_plus,qf_min]),
                    VBox([k_dec,k_acc,k_res,k_tag,k_asymm]),VBox([dm,dg,gs,r,delta,gamma,beta]),
                    VBox([a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t,omega,d_omega,eff,d_eff,a_prod,a_det]),
                    VBox([xmin_coeffs,xmax_coeffs,y_cosh_coeffs,y_sinh_coeffs,y_cos_coeffs,y_sin_coeffs,name_coeffs,save_coeffs])])
wbox_coeffs = HBox([VBox([qt_plus,qt_nul,qt_min,qf_plus,qf_min]),
                    VBox([k_dec,k_acc,k_res,k_tag,k_asymm]),
                    VBox([xmin_coeffs,xmax_coeffs,y_cosh_coeffs,y_sinh_coeffs,y_cos_coeffs,y_sin_coeffs]),
                    VBox([name_coeffs]),
                    VBox([save_coeffs])])
coeffs_pars = interactive_output(plot_coeffs,{
                'qt_plus':qt_plus,'qt_nul':qt_nul,'qt_min':qt_min,'qf_plus':qf_plus,'qf_min':qf_min,
                'dm':dm,'dg':dg,'gs':gs,'r':r,'delta':delta,'gamma':gamma,'beta':beta,
                'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,
                'sigma_t':sigma_t,'omega_tag':omega,'d_omega_tag':d_omega,'eff_tag':eff,'d_eff_tag':d_eff,
                'a_prod_asym':a_prod,'a_det_asym':a_det,
                'xmin':xmin_coeffs,'xmax':xmax_coeffs,'y_cosh':y_cosh_coeffs,'y_sinh':y_sinh_coeffs,'y_cos':y_cos_coeffs,'y_sin':y_sin_coeffs,
                'name':name_coeffs,'save':save_coeffs,
                'k_dec':k_dec,
                'k_acc':k_acc,
                'k_res':k_res,
                'k_tag':k_tag,
                'k_asymm':k_asymm})
