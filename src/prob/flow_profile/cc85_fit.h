#ifndef _FLOW_PROFILE_CC85_FIT_H_
#define _FLOW_PROFILE_CC85_FIT_H_

void flow_profile_init();

Real flow_profile_velocity_norm(Real r);

Real flow_profile_velocity_x(Real x1, Real x2, Real x3);
Real flow_profile_velocity_y(Real x1, Real x2, Real x3);
Real flow_profile_velocity_z(Real x1, Real x2, Real x3);
Real flow_profile_density(Real x1, Real x2, Real x3);
Real flow_profile_pressure(Real x1, Real x2, Real x3);


#endif /* _FLOW_PROFILE_CC85_FIT_H_ */
