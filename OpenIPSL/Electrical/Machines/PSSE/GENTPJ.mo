within OpenIPSL.Electrical.Machines.PSSE;
model GENTPJ
  "WECC Type J GENERATOR: ROUND ROTOR WITH SATURATION ON BOTH AXES."
  extends OpenIPSL.Icons.VerifiedModel;
  //Import of dependencies
  import Complex;
  import Modelica.ComplexMath.arg;
  import Modelica.ComplexMath.real;
  import Modelica.ComplexMath.imag;
  import Modelica.ComplexMath.abs;
  import Modelica.ComplexMath.conj;
  import Modelica.ComplexMath.fromPolar;
  import Modelica.ComplexMath.j;
  import OpenIPSL.NonElectrical.Functions.SE;

  extends OpenIPSL.Electrical.Machines.PSSE.BaseClasses.baseMachine(
    XADIFD(start = efd0), 
    delta(start = 0, fixed = true), 
    id(start = 0), 
    iq(start = 0), 
    Te(start = 0), 
    ud(start = 0), 
    uq(start = 0));

  //Machine parameters
  parameter Types.PerUnit Xpq "q-axis transient reactance "
    annotation (Dialog(group="Machine parameters"));
  parameter Types.Time Tpq0
    "q-axis transient open-circuit time constant"
    annotation (Dialog(group="Machine parameters"));
  parameter Types.PerUnit Kis "Current multiplier for saturation calculation";

  Types.PerUnit Epd(start = 0) "d-axis voltage behind transient reactance ";
  Types.PerUnit Epq(start = 0) "q-axis voltage behind transient reactance ";
  Types.PerUnit Eq1(start = 0);
  Types.PerUnit Eq2(start = 0);
  Types.PerUnit Ed1(start = 0);
  Types.PerUnit Ed2(start = 0);
  Types.PerUnit Xppdsat;
  Types.PerUnit Xppqsat;
  Types.PerUnit dsat;
  Types.PerUnit qsat;
  Types.PerUnit PSId;
  Types.PerUnit PSIq;
  Types.PerUnit PSIppd(start = 0);
  Types.PerUnit PSIppq(start = 0);
  Types.PerUnit PSIpp "Sub-transient flux linkage in stator reference frame";
  Types.PerUnit XadIfd(start = 0);

protected
  parameter Complex VT = fromPolar(v_0, angle_0) "Complex terminal voltage";
  parameter Complex S = Complex(p0, q0) "Complex power on machine base";
  parameter Types.PerUnit vr0 = v_0 * cos(angle_0);
  parameter Types.PerUnit vi0 = v_0 * sin(angle_0);
  parameter Types.PerUnit ir0 = -CoB * (p0 * vr0 + q0 * vi0) / (vr0 ^ 2 + vi0 ^ 2);
  parameter Types.PerUnit ii0 = -CoB * (p0 * vi0 - q0 * vr0) / (vr0 ^ 2 + vi0 ^ 2);
  parameter Types.PerUnit pm0 = real(S);
  parameter Types.PerUnit efd0 = 1.0;
  parameter Types.PerUnit Epq0 = 1.0;
  parameter Types.PerUnit Epd0 = 1.0;
  parameter Types.PerUnit Eq10 = 1.0;
  parameter Types.PerUnit Ed10 = 1.0;
  parameter Types.PerUnit Eq20 = 1.0;
  parameter Types.PerUnit Ed20 = 1.0;
  parameter Types.PerUnit Xppdsat0 = 1.0;
  parameter Types.PerUnit Xppqsat0 = 1.0;
  // Constants
  parameter Real CoB = M_b / S_b "Constant to change from system base to machine base";

initial equation
  der(Epd) = 0;
  der(Epq) = 0;
  der(PSIppd) = 0;
  der(PSIppq) = 0;

equation
  //Interfacing outputs with the internal variable
  XADIFD = XadIfd;
  ISORCE = XadIfd;
  EFD0 = efd0;
  PMECH0 = pm0;

  // Differential equations:
  der(Epq) = (1/Tpd0)*(EFD - XadIfd);
  der(Epd) = (1/Tpq0)*(-1)*qsat*Ed1;
  der(PSIppd) = -(dsat)*((Xpd-Xppd)/(Xd-Xppd))*(Eq2/Tppd0);
  der(PSIppq) = (qsat)*((Xpq-Xppq)/(Xq-Xppq))*(Ed2/Tppq0);
  Te = PSId*iq - PSIq*id;

  // Unsaturated air-gap flux:
  PSIpp = sqrt((uq+iq*R_a+id*Xl)^2 + (ud+id*R_a-iq*Xl)^2);
  // Saturation on d-axis:
  dsat=1+SE((PSIpp + Kis * sqrt(id^2 + iq^2)), S10, S12, 1, 1.2);
  // Saturation on q-axis:
  qsat=1+(Xq/Xd)*SE((PSIpp+Kis*sqrt(id^2+iq^2)), S10, S12, 1, 1.2);

  // Auxiliary Equations:
  Eq1= ((-1)*PSIppd*(Xd-Xpd) + Epq*(Xd-Xppd))/(Xpd-Xppd);
  Ed1= (PSIppq*(Xq-Xpq)+Epd*(Xq-Xppq))/(Xpq-Xppq);
  Eq2= (PSIppd-Epq+id*((Xpd-Xppd)/dsat))*((Xd-Xppd)/(Xpd-Xppd));
  Ed2=-(Epd+PSIppq)*((Xq-Xppq)/(Xpq-Xppq))-iq*((Xq-Xppq)/qsat);

  // Field Current:
  XadIfd = dsat*Eq1;
  // Flux and saturated inductances:
  Xppdsat=((Xppd-Xl)/dsat)+Xl;
  Xppqsat=((Xppq-Xl)/qsat)+Xl;
  PSId=PSIppd - Xppdsat*id;
  PSIq=PSIppq - Xppqsat*iq;
  // Terminal voltage:
  ud = (-PSIq) - R_a*id;
  uq = PSId - R_a*iq;
  annotation (
    Icon( graphics={Text(
          extent={{-54,24},{54,-26}},
          lineColor={0,0,255},
          textString="GENTPJ")}),
    Documentation(info="<html><p>Solid rotor generator with saturation on both axes. The saturation in this model is not only function of the air-gap flux, but also of armature current magnitude. This effect is included via parameter <code>Kis</code>.</p>
    <p> If <code>Kis</code> is set to zero, then the model will behave like the WECC Type F Generator, that is, GENTPF.</p></html>",
    revisions="<html>
<table cellspacing=\"1\" cellpadding=\"1\" border=\"1\">
<tr>
<td><p>Reference</p></td>
<td><p>PSS&reg;E and PowerWorld Manuals</p></td>
</tr>
<tr>
<td><p>Last update</p></td>
<td><p>2025-05-08 Glen Halley</p></td>
</tr>
<tr>
<td><p>Author</p></td>
<td><p>Md. Shamimul Islam and Marcelo de Castro, Rensselaer Polytechnic Institute</p></td>
</tr>
<tr>
<td><p>Contact</p></td>
<td><p>see <a href=\"modelica://OpenIPSL.UsersGuide.Contact\">UsersGuide.Contact</a></p></td>
</tr>
</table>
</html>"));
end GENTPJ;
