within ;

package ADAN_main

  package Components

    package Auxiliary

      package AcausalConnector

        model Pq_terminator_p
          "creates a P type according to Soroushs definition, therefore requires pressure (u) as an input"
          replaceable Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a port_a
            annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(extent={{90,-10},{110,10}})));
          input Physiolibrary.Types.Pressure u;
          Physiolibrary.Types.VolumeFlowRate v = port_a.q;
        equation
          u = port_a.pressure;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                  Rectangle(
                  extent={{-100,100},{-20,-100}},
                  lineThickness=0.5,
                  pattern=LinePattern.None,
                  lineColor={0,0,0},
                  fillColor={244,125,35},
                  fillPattern=FillPattern.Solid), Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={28,108,200},
                  lineThickness=0.5),
                Text(
                  extent={{-86,-42},{-30,40}},
                  lineColor={0,0,0},
                  pattern=LinePattern.None,
                  lineThickness=0.5,
                  fillColor={244,125,35},
                  fillPattern=FillPattern.None,
                  textString="P
type"),         Text(
                  extent={{-20,-100},{100,-20}},
                  lineColor={0,0,0},
                  pattern=LinePattern.None,
                  lineThickness=0.5,
                  fillColor={244,125,35},
                  fillPattern=FillPattern.None,
                  textString="%name")}),                                 Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        end Pq_terminator_p;

        model Pq_terminator_v
          "creates a V type according to Soroushs definition, therefore requires flow (v) as an input"
          replaceable Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a port_a
            annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(extent={{90,-10},{110,10}})));
          Physiolibrary.Types.Pressure u = port_a.pressure;
          input Physiolibrary.Types.VolumeFlowRate v;
        equation
          v = port_a.q;
          annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                  Rectangle(
                  extent={{-100,100},{-20,-100}},
                  lineThickness=0.5,
                  pattern=LinePattern.None,
                  lineColor={0,0,0},
                  fillColor={244,125,35},
                  fillPattern=FillPattern.Solid), Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={28,108,200},
                  lineThickness=0.5),
                Text(
                  extent={{-86,-42},{-30,40}},
                  lineColor={0,0,0},
                  pattern=LinePattern.None,
                  lineThickness=0.5,
                  fillColor={244,125,35},
                  fillPattern=FillPattern.None,
                  textString="V
type"),         Text(
                  extent={{-20,-100},{100,-20}},
                  lineColor={0,0,0},
                  pattern=LinePattern.None,
                  lineThickness=0.5,
                  fillColor={244,125,35},
                  fillPattern=FillPattern.None,
                  textString="%name")}),                                 Diagram(
                coordinateSystem(preserveAspectRatio=false)));
        end Pq_terminator_v;
      end AcausalConnector;

      partial model HeartBase
        extends Physiolibrary.Icons.Heart;


        parameter Boolean UseFrequencyInput = false annotation(choices(checkBox=true));
        parameter Boolean UseThoracicPressureInput = false annotation(choices(checkBox=true));
        parameter Physiolibrary.Types.Frequency HR = 1 "Heart rate, when not specified externally" annotation(Dialog(enabled = not UseFrequencyInput));

        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a sv
          annotation (Placement(transformation(extent={{-110,90},{-90,110}})));
        Physiolibrary.Types.RealIO.FrequencyInput frequency_input if UseFrequencyInput annotation (Placement(
              transformation(extent={{-126,-20},{-86,20}}), iconTransformation(extent={{-120,
                  -20},{-80,20}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a pv
          annotation (Placement(transformation(extent={{-110,-110},{-90,-90}})));
        Physiolibrary.Types.RealIO.PressureInput thoracic_pressure_input if UseThoracicPressureInput  annotation (Placement(
              transformation(extent={{-28,-120},{12,-80}}), iconTransformation(extent={{-20,
                  -120},{20,-80}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b pa
          annotation (Placement(transformation(extent={{90,-110},{110,-90}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b sa
          annotation (Placement(transformation(extent={{90,90},{110,110}})));

        Physiolibrary.Types.Constants.FrequencyConst HR0(k(displayUnit="1/min")=
               HR) if
             not UseFrequencyInput
          annotation (Placement(transformation(extent={{-84,-4},{-76,4}})));
        Physiolibrary.Types.Constants.PressureConst P0(k=0) if
             not UseThoracicPressureInput
          annotation (Placement(transformation(extent={{4,-4},{-4,4}},
              rotation=270,
              origin={0,-80})));
      end HeartBase;

      partial model PulmonaryBase
        extends Physiolibrary.Icons.Lungs;

        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a port_a
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
              iconTransformation(extent={{-110,-10},{-90,10}})));
       Physiolibrary.Types.RealIO.PressureInput thoracic_pressure = _thoracic_pressure if UseThoracic_PressureInput annotation (Placement(
              transformation(extent={{-28,-120},{12,-80}}), iconTransformation(extent={{-20,
                  -120},{20,-80}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b port_b
          annotation (Placement(transformation(extent={{90,-10},{110,10}}),
              iconTransformation(extent={{90,-10},{110,10}})));
        parameter Boolean UseThoracic_PressureInput=true    annotation(choices(checkBox=true));
        Physiolibrary.Types.Pressure _thoracic_pressure;
        Physiolibrary.Types.Volume volume;
      equation
        if not UseThoracic_PressureInput then
          _thoracic_pressure = 0;
        end if;

      end PulmonaryBase;
    end Auxiliary;

    model ConditionalConnection
      extends Modelica.Blocks.Interfaces.SISO;
      Modelica.Blocks.Sources.Constant const1(k=disconnectedValue)
        annotation (Placement(transformation(extent={{-60,-24},{-44,-10}})));
      Modelica.Blocks.Logical.Switch switch1
        annotation (Placement(transformation(extent={{-2,-16},{18,4}})));
      Modelica.Blocks.Sources.BooleanExpression useClosedLoopHR(y=not disconnected)
        annotation (Placement(transformation(extent={{-42,-20},{-22,0}})));
    /*  
  Modelica.Blocks.Interfaces.RealInput u1 annotation (Placement(transformation(
          rotation=0, extent={{-70,-8},{-50,12}}), iconTransformation(extent={{-60,-10},
            {-40,10}})));
  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(
          rotation=0, extent={{52,-10},{72,10}}), iconTransformation(extent={{50,
            -10},{70,10}})));
*/
    parameter Boolean disconnected=false   "When true, the parameter value is used as fixed output" annotation(choices(checkBox=true));
    parameter Real disconnectedValue = 0 "output for disconnected = true. Use SI units!" annotation(Dialog(enable = disconnected));
      MyDelay                              myDelay(   delayTime=delayTime)
        annotation (Placement(transformation(extent={{26,-16},{46,4}})));
      parameter Modelica.SIunits.Time delayTime=0
        "Delay time of output with respect to input signal";
      Modelica.Blocks.Nonlinear.Limiter limiter(uMax=uMax, uMin=uMin)
        annotation (Placement(transformation(extent={{56,-16},{76,4}})));
      parameter Real uMax=Modelica.Constants.inf "Upper limits of input signals";
      parameter Real uMin=-uMax "Lower limits of input signals";
    protected
      Boolean showLimiterLine = not disconnected and (uMax <> Modelica.Constants.inf or uMax <> -uMin);
    equation



      connect(switch1.u3, const1.y) annotation (Line(points={{-4,-14},{-20,-14},{-20,
              -17},{-43.2,-17}}, color={0,0,127}));
      connect(useClosedLoopHR.y,switch1. u2) annotation (Line(points={{-21,-10},{-12,
              -10},{-12,-6},{-4,-6}},         color={255,0,255}));
      connect(y, y) annotation (Line(points={{110,0},{110,0}},
                                                             color={0,0,127}));
      connect(switch1.y, myDelay.u)
        annotation (Line(points={{19,-6},{24,-6}}, color={0,0,127}));
      connect(u, switch1.u1) annotation (Line(points={{-120,0},{-62,0},{-62,2},
              {-4,2}}, color={0,0,127}));
      connect(myDelay.y, limiter.u)
        annotation (Line(points={{47,-6},{54,-6}}, color={0,0,127}));
      connect(limiter.y, y) annotation (Line(points={{77,-6},{90,-6},{90,0},{110,0}},
            color={0,0,127}));
      annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}})),
                                                                         Icon(
            coordinateSystem(extent={{-100,-80},{100,100}}),
                                                          graphics={
            Rectangle(
              extent={{-60,20},{60,-40}},
              pattern=LinePattern.None,
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}),
            Line(points={{-40,0},{52,0}}, color={28,108,200}),
            Rectangle(
              extent={{-10,10},{10,-10}},
              fillColor=DynamicSelect({0,140,72}, if not disconnected then {0,140,72} else {238,46,47}),
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Text(
              visible = DynamicSelect(true, if disconnected then true else false),
              extent={{-60,-34},{60,-14}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.None,
              textString="%disconnectedValue"),
            Line(
              visible = DynamicSelect(true, if showLimiterLine  then true else false),
              points={{-42,-36},{-12,-36},{10,-16},{52,-16}},
              color={238,46,47},
              thickness=0.5)}));
    end ConditionalConnection;

  package Vessel_modules

    package Interfaces

      partial model bg_base
        replaceable Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a port_a annotation (
            Placement(transformation(extent={{-110,-10},{-90,10}}),
              iconTransformation(extent={{-110,-10},{-90,10}})),
              Dialog(tab = "Redeclares"));
        replaceable Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b port_b if not terminator                                                                            annotation (
            Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(
                extent={{90,-10},{110,10}})),
              Dialog(tab = "Redeclares"));
        public
        parameter Boolean terminator = false annotation(choices(checkBox=true));
        parameter Boolean UseInertance = true annotation(choices(checkBox=true));
        parameter Boolean UseOuter_thoracic_pressure=false   annotation(choices(checkBox=true));
        parameter Boolean LimitBackflow=false   "Inserts a one-way valve when true" annotation(choices(checkBox=true),Dialog(group = "Parameters"));

        outer Physiolibrary.Types.Fraction phi "a systemic acitvation fraction, 1 being maximal possible. Normal resting is believed to be 1/4 of the maximum (0.25)";
        parameter Physiolibrary.Types.Fraction phi0 = settings.phi0;

        outer Physiolibrary.Types.Pressure thoracic_pressure;
        parameter Physiolibrary.Types.Fraction thoracic_pressure_ratio=1 annotation (Dialog(enable=UseOuter_thoracic_pressure));

        parameter Boolean UseExercise = false;
        outer Physiolibrary.Types.Fraction Exercise;
        Physiolibrary.Types.Fraction exercise;

        Physiolibrary.Types.Pressure u_in = port_a.pressure;
        Physiolibrary.Types.VolumeFlowRate v_in = port_a.q;

        Physiolibrary.Types.Pressure u_out_valved = port_b.pressure if not terminator;
        Physiolibrary.Types.VolumeFlowRate v_out_valved = port_b.q if not terminator;

        Physiolibrary.Types.Pressure u_out;
        Physiolibrary.Types.VolumeFlowRate v_out;

        Physiolibrary.Types.Volume volume(nominal = 1e-6);
        parameter Physiolibrary.Types.Volume V_min = 0;
      /*
protected 
              Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b port_b_internal if not 
    terminator
    annotation (Placement(transformation(extent={{10,-10},{30,10}}),
        iconTransformation(extent={{30,-10},{50,10}})), Dialog(tab="Redeclares"));
public 
  Physiolibrary.Hydraulic.Components.IdealValve idealValve if LimitBackflow
    annotation (Placement(transformation(extent={{40,-20},{60,0}})));
/*  Auxiliary.OnePort_None onePort_None if not LimitBackflow
    annotation (Placement(transformation(extent={{40,0},{60,20}})));
    */
        Real passableVariable(start=0, final unit="1")
          "Auxiliary variable for actual position on the ideal diode characteristic";
        parameter Physiolibrary.Types.HydraulicResistance R_on(final min=0, displayUnit="l/(mmHg.min)") = 0
          "Forward state-on conductance (open valve resistance)" annotation (Dialog(enable=LimitBackflow));
        parameter Physiolibrary.Types.HydraulicConductance G_off(final min=0, displayUnit="l/(mmHg.min)") = 1.2501026264094e-12
          "Backward state-off conductance (closed valve conductance)" annotation (Dialog(enable=LimitBackflow));
        parameter Physiolibrary.Types.Pressure Pknee(final min=0) = 0
          "Forward threshold pressure" annotation (Dialog(enable=LimitBackflow));
        Boolean open(start = true);
        Physiolibrary.Types.Pressure dp = u_out - u_out_valved;
        outer Settings settings
          annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
        protected
        constant Physiolibrary.Types.Pressure unitPressure=1;
        constant Physiolibrary.Types.VolumeFlowRate unitFlow=1;
      equation
       if UseExercise then
         exercise = Exercise;
       else
         exercise = 0;
       end if;


        v_out + v_out_valved = 0;

        if LimitBackflow then
          open = passableVariable > Modelica.Constants.eps;
          dp = (passableVariable*unitFlow)*(if open then R_on else 1) + Pknee;
          v_out = (passableVariable*unitPressure)*(if open then 1 else G_off) + G_off*Pknee;
        else
          // the valve is nonexistent, therefore permanently open
          open = true;
          dp = 0;
          passableVariable = 0;
        end if;

        assert(volume > 0, "Volume is negative!", AssertionLevel.warning);
      /*
  connect(port_b_internal, idealValve.q_in) annotation (Line(
      points={{20,0},{40,0},{40,-10}},
      color={0,0,0},
      thickness=1));
  connect(port_b, idealValve.q_out) annotation (Line(
      points={{100,0},{96,0},{96,-10},{60,-10}},
      color={0,0,0},
      thickness=1));
/*  connect(port_b_internal, onePort_None.q_in) annotation (Line(
      points={{20,0},{30,0},{30,10},{40,10}},
      color={0,0,0},
      thickness=1));
  connect(port_b, onePort_None.q_out) annotation (Line(
      points={{100,0},{80,0},{80,10},{60,10}},
      color={0,0,0},
      thickness=1));
      */
          annotation (Icon(coordinateSystem(extent={{-100,-20},{100,20}}), graphics={
              Text(
                extent={{-100,-20},{100,0}},
                lineColor={28,108,200},
                textString="%name"),
              Rectangle(extent={{-100,20},{100,-20}}, lineColor={28,108,200}),
                Rectangle(extent={{-100,20},{100,-20}}, lineColor={0,140,72},
                  lineThickness =                                                          1,
                visible = DynamicSelect(false, UseOuter_thoracic_pressure)),
              Text(
                extent={{-100,-40},{100,-20}},
                lineColor={0,140,72},
                textString="%thoracic_pressure_ratio",
                visible = DynamicSelect(true, thoracic_pressure_ratio <> 1)),
              Line(
                points={{20,40},{50,10},{80,4}},
                color={28,108,200},
                smooth=Smooth.Bezier,
                thickness=0.5,
                visible = DynamicSelect(true, LimitBackflow)),
              Line(
                points={{20,-40},{50,-10},{80,-4}},
                color={28,108,200},
                smooth=Smooth.Bezier,
                thickness=0.5,
                visible = DynamicSelect(true, LimitBackflow))}),                 Diagram(coordinateSystem(extent={{-100,
                  -20},{100,20}})));
      end bg_base;

      partial model bg_vessel
        extends bg_base(final terminator=false);
        outer Modelica.SIunits.Angle Tilt;
        outer Physiolibrary.Types.Fraction Exercise;
        Physiolibrary.Types.Pressure P_hs = sin(Tilt)*height*rho*Modelica.Constants.g_n "Hydrostatic pressure";
        Physiolibrary.Types.Pressure u_out_hs = u_out + P_hs "Output pressure including the hydrostatic pressure";

        parameter Boolean UseDistentionOutput = false "Provides relative distention fraction output, otherwise hidden and not calculated" annotation(choices(checkBox=true));
        parameter Boolean CalculateMeans = false "uses additional calculations" annotation(choices(checkBox=true));
      //  parameter Boolean UseNonLinearCompliance = false annotation(Dialog(group = "Parameters"));
        parameter Boolean UseSinAlphaInput = false "Vessel vertical angle, 0 is supine and 1 is facing upwards, -1 downwards respectively"  annotation(choices(checkBox=true),Dialog(group = "Parameters"));
        parameter Boolean UseVasoconstrictionEffect = false annotation(Dialog(group = "Parameters"));

        parameter Real a(unit = "1") = 0.2802 annotation (Dialog(tab = "Defaults", group = "Vessel parameter computations"));
        parameter Real b(unit = "m-1") = -505.3 annotation (Dialog(tab = "Defaults", group = "Vessel parameter computations"));
        parameter Real c(unit = "1") = 0.1324 annotation (Dialog(tab = "Defaults", group = "Vessel parameter computations"));
        parameter Real d(unit = "m-1") = -11.14 annotation (Dialog(tab = "Defaults", group = "Vessel parameter computations"));

        constant Real mu(unit = "J.s.m-3") = 0.004;
        constant Real rho(unit = "J.s2.m-5") = 1050;
        parameter Real sinAlpha = 0 "sin of vessel orientation angle, 0 being supine, 1 being up, -1 aiming down."  annotation (Dialog(tab = "General", group = "Orientations", enabled = not UseSinAlphaInput));
        Modelica.SIunits.Height height = _sinAlpha*l annotation (Dialog(tab = "General", group = "Orientations"));

        //parameter Physiolibrary.Types.Pressure thoracic_pressure = 0;

        parameter Real E(unit = "Pa") = 4e5 "Elasticity"  annotation (Dialog(tab = "General", group = "Vessel properties", enable = not UseNonLinearCompliance));
       parameter Modelica.SIunits.Length l = 1e-2 "Segmant length" annotation (Dialog(tab = "General", group = "Vessel properties"));
       parameter Modelica.SIunits.Radius r = 1e-3 "Vessel radius" annotation (Dialog(tab = "General", group = "Vessel properties"));

        parameter Modelica.SIunits.Thickness h = r*(a*exp(b*r)+c*exp(d*r)) "Thickness" annotation (Dialog(tab = "General", group = "Vessel properties"));


        parameter Physiolibrary.Types.HydraulicInertance I = rho*l/(Modelica.Constants.pi*(r)^2) annotation (Dialog(tab = "General", group = "Calculated parameters"));
        parameter Physiolibrary.Types.HydraulicCompliance C = 2*Modelica.Constants.pi*(r^3) *l/(E*h) annotation (Dialog(tab = "General", group = "Calculated parameters", enable = not UseNonLinearCompliance));
        parameter Physiolibrary.Types.HydraulicResistance R = 8*mu*l/(Modelica.Constants.pi*(r^4)) annotation (Dialog(tab = "General", group = "Calculated parameters"));
        parameter Physiolibrary.Types.HydraulicResistance R_v = 0.01/C "Viscoleasticity of the vessel" annotation (Dialog(tab = "General", group = "Calculated parameters"));

        parameter Physiolibrary.Types.Volume zpv = l*Modelica.Constants.pi*(r^2) "Zero-pressure volume" annotation (Dialog(tab = "General", group = "Calculated parameters"));
        Physiolibrary.Types.RealIO.FractionOutput distentionFraction = sqrt(max(volume, 0))/sqrt(distentionBase) if
          UseDistentionOutput annotation (Placement(transformation(extent={{76,10},{96,
                  30}}), iconTransformation(extent={{-20,-20},{20,20}},
              rotation=90,
              origin={0,40})));
        parameter Physiolibrary.Types.Volume distentionBase = l*Modelica.Constants.pi*(r^2);

        Modelica.Blocks.Interfaces.RealInput sinAlphaInput = _sinAlpha if
          UseSinAlphaInput annotation (Placement(transformation(extent={{-84,10},{-64,
                  30}}), iconTransformation(extent={{-20,-20},{20,20}},
              rotation=90,
              origin={-80,-20})));

        Physiolibrary.Types.HydraulicCompliance compliance "Real compliance after all effects";
        Physiolibrary.Types.Fraction vc_effect = 1 - settings.R_vc*(phi - phi0);

        protected
        Real _sinAlpha;
      equation
        if not UseSinAlphaInput then
          _sinAlpha = sinAlpha;
        end if;

        if UseVasoconstrictionEffect then
          compliance = 2*Modelica.Constants.pi*((r*vc_effect)^3) *l/(E*h)/exp(Exercise*settings.exercise_factor_on_arterial_compliance);
        else
          compliance = C;
        end if;
      //  h = r*(a*exp(b*r)+c*exp(d*r));
      //  I = rho*l/(Modelica.Constants.pi*(r)^2);
      //  C = 2*Modelica.Constants.pi*(r^3) *l/(E*h);
      //  R = 8*mu*l/(Modelica.Constants.pi*(r^4));
      //  R_v = 0.01/C;
        annotation (Icon(graphics={
              Text(
                extent={{-100,0},{100,20}},
                lineColor={28,108,200},
                textString=DynamicSelect("L, d",
                "L = " + String(l*100, significantDigits=2) + "cm, " +
                "D = " + String(r*2*100, significantDigits=2) + "cm")),
                Line(
                points=DynamicSelect({{-60,20},{-20,40}}, {{-60,20},{-60 + cos(asin(sinAlpha))*40, 20 + sinAlpha*40}}),
                color={0,140,72},
                thickness=1,
                visible = DynamicSelect(false, sinAlpha > 0 or sinAlpha < 0),
                arrow={Arrow.None,Arrow.Filled})}));
      end bg_vessel;

      partial model bg_vessel_thoracic
        extends bg_vessel;
        annotation (Icon(graphics={Rectangle(
                extent={{-100,20},{100,-20}},
                lineColor={0,140,72},
                lineThickness=0.5)}));
      end bg_vessel_thoracic;

      partial model systemic_tissue_base
        extends ADAN_main.Components.Vessel_modules.Interfaces.bg_base(UseInertance = false, volume(start = V_n, fixed = false));

        parameter Real I(unit = "J.s2.m-6");
        parameter Real C(unit = "m6.J-1");
        parameter Real Ra(unit="J.s.m-6") "Arteriole resistance";
        Real Rvis(unit="J.s.m-6") "Elastic viscosity using Voigt model of in-series resistance";
        Real I_e(unit = "J.s2.m-6");
        parameter Real Rv(unit="J.s.m-6") "venule resistance";

        parameter Physiolibrary.Types.Volume zpv = 0 "Zero-pressure volume";
        parameter Physiolibrary.Types.Pressure nominal_pressure = settings.tissues_nominal_pressure;
        Physiolibrary.Types.Pressure u_C(start = nominal_pressure, nominal = 1000, fixed = true);

        Physiolibrary.Types.Pressure u(nominal = 1000);

        Physiolibrary.Types.Pressure u_out_hs "Output pressure including the hydrostatic pressure";

        Physiolibrary.Types.HydraulicResistance Ra_phi_inf = Ra*exp((phi-phi0)*settings.Ra_factor)/exp(exercise*settings.exercise_factor) "Arterioles resistance dependent on phi";
        Physiolibrary.Types.HydraulicResistance Ra_phi(start = Ra, fixed = true) "Delayed arterioles resistance dependent on phi";

        Physiolibrary.Types.HydraulicResistance Rv_phi = Rv*exp((phi-phi0)*settings.Rv_factor) "Arterioles resistance dependent on phi";

        parameter Real k = C / (V_max - V_n) "For Pstras non-linear PV characteristics";
        parameter Physiolibrary.Types.Volume V_max = V_n + (V_n - zpv)*settings.tissues_gamma
          " V_n is between zpv and V_max, parameterized by gamma. 1 means its in the center. For Pstras non-linear PV characteristics";
        parameter Physiolibrary.Types.Volume V_n = nominal_pressure*C + zpv
          "nominal volume calculated from the linear relationship. nominal volume for non-linear compliance parametrization";
        parameter Physiolibrary.Types.Volume V_us = V_max - (V_max - V_n)*exp(nominal_pressure * C / (V_max - V_n))
          "nonlinear Un-Stressed volume";

        Real k_phi = C / (V_max - V_n_phi) "For Pstras non-linear PV characteristics";
        Physiolibrary.Types.Volume V_max_phi = V_max - (V_us - V_us_phi) "From Pstras";
        Physiolibrary.Types.Volume V_n_phi = V_n * (1 - settings.tissuesCompliance_PhiEffect*(phi - settings.phi0))/exp(exercise*settings.exercise_factor_on_arterial_compliance) "Linearly dependent on phi";
        Physiolibrary.Types.Volume V_us_phi = V_us * (1 - settings.tissuesCompliance_PhiEffect*(phi - settings.phi0)) "Linearly dependent on phi";

      //   Physiolibrary.Types.Fraction phi_shift=(1 + settings.tissues_compliance_phi_shift
      //       *(phi - phi0)) "Nonlinear compliance unstressed volume affected by phi";


      initial equation
      //  volume = nominal_pressure*C + zpv;
      equation

          der(Ra_phi)*settings.tissues_Ra_tau =  Ra_phi_inf - Ra_phi;

            I_e = I*1e-6;
            Rvis = 0.01/C;

            if UseInertance then
              der(v_in) =(u_in - u - Ra_phi*v_in)/I;
              der(v_out) =(u - u_out_hs - Rv_phi*v_out)/I_e;
            else
              0 =(u_in - u - Ra_phi *v_in);
              0 =(u - u_out_hs - Rv_phi *v_out);
            end if;

            der(volume) = (v_in-v_out);
            if UseOuter_thoracic_pressure then
              u =u_C + Rvis*(v_in - v_out) + thoracic_pressure;
            else
              u =u_C + Rvis*(v_in - v_out);
            end if;

            if not settings.UseNonLinear_TissuesCompliance then
              volume = (u_C) *C + zpv;
            elseif settings.tissues_UseStraighteningReaction2Phi then
              // my expression based on straightening
              u_C = 1/k_phi*log((V_max - V_us)/max(1e-9, V_max - volume)) "equation from Pstras 2017, 10.1093/imammb/dqw008, allegedly taken from Hardy & Collins 1982";

            else
              // from Pstras
      //        u_C = 1/k*log((V_max - V_us)/(V_max - volume)) "equation from Pstras 2017, 10.1093/imammb/dqw008, allegedly taken from Hardy & Collins 1982";
              u_C = 1/k*log((V_max_phi - V_us_phi)/max(1e-9, V_max_phi - volume)) "equation from Pstras 2017, 10.1093/imammb/dqw008, allegedly taken from Hardy & Collins 1982";

            end if;

        annotation (Icon(graphics={Polygon(
                points={{-20,-26},{-2,2},{-18,2},{12,30},{8,30},{24,36},{30,30},{26,30},
                    {6,10},{18,10},{-6,-24},{-20,-26}},
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid,
                pattern=LinePattern.None,
                visible=DynamicSelect(false, UseExercise)),
              Rectangle(
                extent={{-20,20},{20,0}},
                lineThickness=0.5,
                fillColor={244,125,35},
                fillPattern=FillPattern.Solid,
                pattern=LinePattern.None)}), Documentation(info="<html>
<p>When <span style=\"font-family: Courier New;\">UseNonLinearCompliance</span> is true, then this formula is applied (from <span style=\"font-family: Courier New; color: #006400;\">Pstras, Math Med Biol 2017,&nbsp;10.1093/imammb/dqw008):</span></p>
<p><img src=\"modelica://ADAN_main/Resources/tissuePV.png\"/></p>
<p>Otherwise, the linear (green) relation is used instead.</p>
<p><br>Unlike the original Pstras publication, when the VMax is set as Vn + 500ml for veins, resp. Vn + 100 ml for Vena Cava, the Vmax difference here is delta-times the ZPV-Vn:</p>
<p><img src=\"modelica://ADAN_main/Resources/Images/equations/equation-sfHOF1Pk.png\" alt=\"Vmax = Vn + delta*(Vn-ZPV)\"/> </p>
<p>and thus shares the linear parametrization.</p>
</html>"));
      end systemic_tissue_base;

      partial model compliance_base
        "A base model for pressure - volume characteristics"

        import Physiolibrary.Types.*;
        input Volume V "Volume";
        outer Fraction phi;
        parameter Fraction phi0 = settings.phi0;
        Fraction phi_ = if settings.veins_UsePhiEffect then phi else phi0 "internal on-off phi";
      //  parameter Boolean UsePhiEffect = true;
        parameter Modelica.SIunits.Length l;
        parameter Volume V0 "Initial volume";
        parameter Volume V_min "minimal collapsing pressure";
        parameter Pressure p0 "Nominal pressure for initialization";
        Pressure p "Fluid pressure";


        outer Settings settings
          annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end compliance_base;

      model compliance_tensionBased
        "Exponential tension-based compliance with linear tension muscle activation. Idea based on Carlson, Secomb 2015"
        extends compliance_base(
          V(start=V0, fixed=true),
          V0=l*Modelica.Constants.pi*(r_n)^2,
          V_min=Modelica.Constants.pi*(wall_L_min/Modelica.Constants.pi/2)^2*l);
        parameter Modelica.SIunits.Radius r_n "Nominal vessel inner radius";
        Modelica.SIunits.Diameter D(start=2*r_n) "Vessel actual (inner) diameter";

        Physiolibrary.Types.Fraction ll0( start = 1)=wall_L/wall_L0;
        Modelica.SIunits.Length wall_L=Modelica.Constants.pi*D
          "Circumferential wall length";
        parameter Modelica.SIunits.Length wall_L0=Modelica.Constants.pi*r_n*2
          "Circumferential wall length at nominal";
        parameter Modelica.SIunits.Length wall_L_min=wall_L0*gamma "Minimal circumferential length";
        parameter Real gamma = settings.veins_gamma "fraction of nominal diameter to minimal zero-pressure diameter";
        parameter Real alpha = settings.veins_alpha "how many times the tension is larger for maximal activation from resting activation at nominal diameter";

        Physiolibrary.Types.Fraction A_inf = phi_ "Instant muscle activation coefficient";
        Physiolibrary.Types.Fraction A "Delayed muscle activation coefficient";
        parameter Physiolibrary.Types.Time tau=5   "Time constant of the smooth muscle activation";

        Real T_total = T_pass + T_active "Total vessel wall tension";
        Real T_pass = C_pass*(exp(T_pass_exp) - 1) "Passive part of the waqll tension";
        Real T_pass_exp= (wall_L - wall_L_min) /wall_L_min;
        Real T_active = A*T_act_max "Active part of the wall tension";
        Real T_act_max= C_act*(wall_L - wall_L_min)/wall_L_min "Tension at maximal muscle activation";

      // Calculation of primary parameters
        Real T_pass_base = C_pass*(exp((wall_L0 - wall_L_min) /wall_L_min) - 1) "passive tension at nominal l = l0";
        Real T_nominal = T_pass_base + C_act*phi0 "Total tension at nominal (l = l0) and resting activation";
        Real T_nominal_max=T_pass_base + C_act
          "Total tension at nominal and maximal activation";


        Real C_pass(start = 1) "Primary parameter for passive tension";
        Real C_act(start = 1) "Primary parameter for active tension";
        Real C_act_max=p0*r_n*4;
        Real C_act_guess=p0*r_n*(gamma*(alpha - 1)/(1 + phi_*(alpha - 1)))
          "Analytically expressed equation";
        Real C_pass_guess=(p0*r_n - C_act/4)/(exp((1 - gamma)/gamma) - 1);


      initial equation
        A = A_inf;
      equation
        assert(C_pass > 0, "C_pass is negative - change alpha or gamma parameters");

        V = l * Modelica.Constants.pi * (D/2)^2;
        der(A)*tau = A_inf - A;

        // identify the c_act
        alpha =T_nominal_max / T_pass_base;
        // identify the c_pass - the vessel must exhibit the given p0 at nominal diameter
        p0*r_n = T_nominal;


        // governing equation - tension to pressure
        p * D/2 = T_total;

      end compliance_tensionBased;

      model compliance_dataFit1
        "Fit of venous PV data from Moreno 1970 by Ben (ebrandal@umich.edu)"
        extends compliance_base(
          V0 = Modelica.Constants.pi*r_n^2*l,
          V_min = Modelica.Constants.pi*r_0^2*l);
        type Tension =   Real (final quantity = "Tension", final unit = "N/m");
        type Radius =   Modelica.SIunits.Radius(final nominal =   1e-6);
        type CircumferentialLength =   Modelica.SIunits.Length(final nominal =   1e-3);
        // GENERAL INPUT PARAMETERS
      //  parameter Boolean useViscoElasticDelay = false;
      //   parameter Physiolibrary.Types.Fraction gamma =   0.5
      //     "Fraction of minimal collapsing diameter to nominal diameter";

        // TENSIONS
        Tension T = T_p + T_a;
        Tension T_p = a*f_L + b*g_L "Passive vessel wall tension";
        Tension T_a_max = c*h_L "Maximal active vessel wall tension";
        Tension T_a = A*T_a_max "Actual active vessel wall tension";
        Real f_L = L*(L - L_0)/L_0 "Tension function on circumferential wall length";
        Real g_L = L_0*(exp(d*(L - L_0)/L_0) - 1)
          "Tension function on circumferential wall length";
        Real h_L = (L - L_0)*exp(-e*((L - L_0)/L_0)^2)
          "Tension function on circumferential wall length";

        // TENSION PARAMETERS
        constant Real mmHg2Pa = 133.32;
        parameter Tension a = mmHg2Pa*1.2 "identified by Matlab's cftool and fixed";
        parameter Tension b = mmHg2Pa*6e-5 "identified by Matlab's cftool and fixed";
        parameter Tension c = mmHg2Pa*26.25 "identified by Matlab's cftool and fixed";
        parameter Real d = 11.5 "identified by Matlab's cftool and fixed";
        parameter Real e = 1.5;

        // DIAMETERS AND LENGTHS
        Radius r(start = r_n, fixed = false) "Actual vessel radius";
        parameter Radius r_n "nominal vessel radius";
        parameter Radius r_0 = sqrt(V_0/(l*Modelica.Constants.pi)) "Zero-pressure vessel radius";
        CircumferentialLength L = r*(2*Modelica.Constants.pi)
          "Actual circumferential wall length";
        parameter CircumferentialLength L_n = r_n*(2*Modelica.Constants.pi)
          "nominal circumferential wall length";
        parameter CircumferentialLength L_0 = r_0*(2*Modelica.Constants.pi)
          "Zero-pressure circumferential wall length";
        parameter Physiolibrary.Types.Volume V_0 =   V_n*settings.veins_gamma "Zero-pressure volume";
        parameter Physiolibrary.Types.Volume V_n =   V0 "nominal aka initial volume";

        // ACTIVE REGULATION
        Physiolibrary.Types.Fraction A(start=A_n) "Activation fraction";
        parameter Physiolibrary.Types.Fraction A_n=phi0 "nominal activation fraction";
        parameter Modelica.SIunits.Time tau = settings.veins_activation_tau;


        // INITIALIZATION - calculations of the nominal tensions and pressures
        Tension T_n = a*f_L_n + b*g_L_n + A_n*c*h_L_n;
        Real f_L_n = L_n*(L_n - L_0)/L_0 "Tension function on circumferential wall length";
        Real g_L_n = L_0*(exp(d*(L_n - L_0)/L_0) - 1)
          "Tension function on circumferential wall length";
        Real h_L_n = (L_n - L_0)*exp(-e*((L_n - L_0)/L_0)^2)
          "Tension function on circumferential wall length";

        Physiolibrary.Types.Pressure p_n;

        // ASSUMPTIONS used for parameter identification
        //  parameter Tension T_n  =  p0 * r_n "Tension at nominal pressure of p0";
        //  parameter Physiolibrary.Types.Pressure P_dm  =  30*133.32 "Data point - maximal pressure";
        //  parameter Physiolibrary.Types.Volume V_dm  =  4*V0 "Volume data point, corresponding to P_dm";
        //   parameter Modelica.SIunits.Length L_dm  =  2*Modelica.Constants.pi*r_dm;
        //   parameter Modelica.SIunits.Radius r_dm  =  sqrt(V_dm/Modelica.Constants.pi/l);
        //   parameter Tension T_dm  =  P_dm * r_dm;
        //   parameter Physiolibrary.Types.Fraction alpha  =  2.5;

        //   function f
        //     input Modelica.SIunits.Length L_i;
        //     input Modelica.SIunits.Length L0;
        //     output Tension T;
        //   algorithm
        //     T : =  L_i*(L_i - L0)/ L0^2;
        //   end f;
        //
        //   function g
        //     input Modelica.SIunits.Length L_i;
        //     input Modelica.SIunits.Length L0;
        //     output Tension T;
        //   algorithm
        //     T : =  exp(c*(L_i-L0)/L0) - 1;
        //   end g;
        //
        //   function h
        //     input Modelica.SIunits.Length L_i;
        //     input Modelica.SIunits.Length L0;
        //     output Tension T;
        //   algorithm
        //     T : =  (L_i-L0) / L0;
        //   end h;

        //   Real a = (T_n/(1 + A_nominal*(alpha - 1)) - b*g(L_n, L_0))/f(L_n, L_0);
        //   Real b = (T_dm - T_n/(1 + A_nominal*(alpha - 1))*f(L_dm, L_0)/f(L_n, L_0))/(g(
        //       L_dm, L_0) - g(L_n, L_0)*f(L_dm, L_0)/f(L_n, L_0));
        //   Real d = T_n*(alpha - 1)/(h(L_n, L_0)*(1 + A_nominal*(alpha - 1)));
        //  Real helper  =  T_n/(1+ A_nominal*(alpha -1));
        //  Real helper2  =  f(L_dm)/f(L_n);

      equation


        der(A)*tau  =  phi_ - A;

        // VOLUME equation
        Modelica.Constants.pi*r^2*l  =  V;

        // pressure-tension equation
        p*r  =  T;

        // nominal pressure-tension equatino
        p_n*r_n = T_n;

      end compliance_dataFit1;
    end Interfaces;

    model vv_type_thoracic
      extends Interfaces.bg_vessel_thoracic;
    //  input Physiolibrary.Types.Pressure thoracic_pressure;

    //  Real u(unit = "Pa");
    //  Real u_C(unit = "Pa", start = 0.0);
      Real v(unit = "m3.s-1", start = 0.0);
      Real u_d(unit = "Pa");
      Real u_C_d(unit = "Pa", start = 0.0);
    equation
        volume = u_in*C/2 + u_C_d*C/2 + zpv;
        //u_in = u_C;
       // u = u_out;

          der(v) = (u_out-u_d-R*v)/I;
          der(u_in - thoracic_pressure) = (v_in-v)/(C/2);
          der(u_C_d) = (v-v_out)/(C/2);
          u_out = u_in+2*R_v*(v_in-v);
          u_d = u_C_d+2*R_v*(v-v_out);

      annotation (Icon(graphics={
            Line(
              points={{-100,0},{-60,0}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Open}),
            Line(
              points={{60,0},{100,0}},
              color={28,108,200},
              arrow={Arrow.None,Arrow.Open})}));
    end vv_type_thoracic;

    model pv_type
        extends ADAN_main.Components.Vessel_modules.Interfaces.bg_vessel(
            UseVasoconstrictionEffect=settings.arteries_UseVasoconstrictionEffect,
            UseExercise = true);

      Real u_C(unit = "Pa", start = 10000.0, fixed = true);



    equation
      volume = (u_C)  *compliance + zpv "Lim 2013";

      if UseInertance then
        der(v_in) = (u_in-u_out_hs-R*v_in)/I;
      else
        0 = u_in-u_out_hs-R*v_in;
      end if;

      der(volume) = (v_in-v_out);
      if UseOuter_thoracic_pressure then
          u_out_hs = u_C+R_v*(v_in-v_out) + thoracic_pressure;
      else
        u_out_hs = u_C+R_v*(v_in-v_out);
      end if;

        annotation (Icon(graphics={Line(
                points={{-80,0},{80,0}},
                color={238,46,47},
                arrow={Arrow.None,Arrow.Filled},
                thickness=0.5)}));
    end pv_type;

    model pv_type_thoracic
      extends ADAN_main.Components.Vessel_modules.pv_type(final
            UseOuter_thoracic_pressure=true);
    end pv_type_thoracic;

    model vp_type

      extends Interfaces.bg_vessel(
          UseOuter_thoracic_pressure=false,
      UseInertance = false,
      zpv = l*Modelica.Constants.pi*((r*venous_diameter_correction)^2),
      R = 8*mu*l/(Modelica.Constants.pi*((r*venous_diameter_correction)^4)),
      I = rho*l/(Modelica.Constants.pi*(r*venous_diameter_correction)^2),
      V_min = compliant_vessel.V_min,
      volume(start = compliant_vessel.V0, fixed = true));

    //   outer Modelica.SIunits.Angle Tilt;
      parameter Physiolibrary.Types.Fraction venous_diameter_correction = settings.venous_diameter_correction;

      parameter Physiolibrary.Types.Pressure p0= settings.tissues_nominal_venules_pressure "nominal venous pressure";
    //  parameter Boolean ignoreViscosityResistance = true;
    //  parameter Boolean limitExternalPressure =  true;


    //   Physiolibrary.Types.Pressure P_hs = sin(Tilt)*height*rho*Modelica.Constants.g_n "Hydrostatic pressure";
    //   Physiolibrary.Types.Pressure u_out_hs = u_out + P_hs "Output pressure including the hydrostatic pressure";

      Physiolibrary.Types.Pressure external_pressure;
      Physiolibrary.Types.VolumeFlowRate netFlow = v_in-v_out;

      Physiolibrary.Types.Pressure p = compliant_vessel.p;


        replaceable Interfaces.compliance_tensionBased compliant_vessel(
          l=l,
          r_n=r) constrainedby Interfaces.compliance_base(
          l=l,
          p0=p0,
          V=volume)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation

      if E == 0 or not UseInertance then
          (u_in-u_out_hs)/R = v_out;
      else
        der(v_out) = (u_in-u_out_hs-R*v_out)/I;
      end if;

      der(volume) = netFlow;


      if UseOuter_thoracic_pressure and settings.veins_limitExternalPressure then
        external_pressure =
          thoracic_pressure_ratio*thoracic_pressure
          *max(min(compliant_vessel.V /compliant_vessel.V_min, 1), 0);
      elseif UseOuter_thoracic_pressure then
        external_pressure = thoracic_pressure_ratio*thoracic_pressure;
      else
        external_pressure = 0;
      end if;

      if not settings.veins_ignoreViscosityResistance then
          u_in = p + R_v*(v_in - v_out) + external_pressure;
        else
          u_in = p + external_pressure;
      end if;

    end vp_type;

    model vp_type_tension_based "Extension of datafit"
      extends Interfaces.bg_vessel(
          UseOuter_thoracic_pressure=false,
      UseInertance = false,
      zpv = l*Modelica.Constants.pi*((r*venous_diameter_correction)^2),
      R = 8*mu*l/(Modelica.Constants.pi*((r*venous_diameter_correction)^4)),
      I = rho*l/(Modelica.Constants.pi*(r*venous_diameter_correction)^2),
      V_min = compliant_vessel.V_min,
      volume(start = compliant_vessel.V0, fixed = true));

    //   outer Modelica.SIunits.Angle Tilt;
      parameter Physiolibrary.Types.Fraction venous_diameter_correction = settings.venous_diameter_correction;

      parameter Physiolibrary.Types.Pressure p0= settings.tissues_nominal_venules_pressure "nominal venous pressure";
    //  parameter Boolean ignoreViscosityResistance = true;
    //  parameter Boolean limitExternalPressure =  true;


    //   Physiolibrary.Types.Pressure P_hs = sin(Tilt)*height*rho*Modelica.Constants.g_n "Hydrostatic pressure";
    //   Physiolibrary.Types.Pressure u_out_hs = u_out + P_hs "Output pressure including the hydrostatic pressure";

      Physiolibrary.Types.Pressure external_pressure;
      Physiolibrary.Types.VolumeFlowRate netFlow = v_in-v_out;

      Physiolibrary.Types.Pressure p = compliant_vessel.p;


        Interfaces.compliance_dataFit1 compliant_vessel(
          l=l,
          r_n=r,
          p0=p0,
          V=volume)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation

      if E == 0 or not UseInertance then
          (u_in-u_out_hs)/R = v_out;
      else
        der(v_out) = (u_in-u_out_hs-R*v_out)/I;
      end if;

      der(volume) = netFlow;


      if UseOuter_thoracic_pressure and settings.veins_limitExternalPressure then
        external_pressure =
          thoracic_pressure_ratio*thoracic_pressure
          *max(min(compliant_vessel.V /compliant_vessel.V_min, 1), 0);
      elseif UseOuter_thoracic_pressure then
        external_pressure = thoracic_pressure_ratio*thoracic_pressure;
      else
        external_pressure = 0;
      end if;

      if not settings.veins_ignoreViscosityResistance then
          u_in = p + R_v*(v_in - v_out) + external_pressure;
        else
          u_in = p + external_pressure;
      end if;
    end vp_type_tension_based;

    model systemic_tissue
      extends
          ADAN_main.Components.Vessel_modules.Interfaces.systemic_tissue_base;
    equation

        u_out_hs = u_out;

      annotation (Icon(graphics={
            Rectangle(
              extent={{-20,20},{20,0}},
              lineThickness=0.5,
              fillColor={244,125,35},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None)}));
    end systemic_tissue;
  end Vessel_modules;

    package AdanVenousRed

      package _b580e "Commit to 4ac repository from May 2019"

      package Parameters_Venous_cellml

        model Parameters_Pulmonary
          parameter Real C_pas(unit = "m6.J-1") = 0.00135e-6;
          parameter Real C_pat(unit = "m6.J-1") = 0.0285e-6;
          parameter Real C_par(unit = "m6.J-1") = 0.0232e-6;
          parameter Real C_pcp(unit = "m6.J-1") = 0.0684e-6;
          parameter Real C_pvn(unit = "m6.J-1") = 0.15376e-6;
          parameter Real C_pvc(unit = "m6.J-1") = 0.01125e-6;
          parameter Real R_pas(unit = "J.s.m-6") = 0.26664e+6;
          parameter Real R_pat(unit = "J.s.m-6") = 1.3332e+6;
          parameter Real R_par(unit = "J.s.m-6") = 6.666e+6;
          parameter Real R_pcp(unit = "J.s.m-6") = 33.33e+6;
          parameter Real R_psh(unit = "J.s.m-6") = 701.11e+6;
          parameter Real R_pvn(unit = "J.s.m-6") = 0.8e+6;
          parameter Real I_pas(unit = "J.s2.m-6") = 0.00693e+6;
          parameter Real I_pat(unit = "J.s2.m-6") = 0.226644e+6;
          parameter Real I_par(unit = "J.s2.m-6") = 1e-12;
          parameter Real I_pcp(unit = "J.s2.m-6") = 1e-12;
          parameter Real I_pvn(unit = "J.s2.m-6") = 1e-12;
          parameter Real I_psh(unit = "J.s2.m-6") = 1e-12;
        equation

        end Parameters_Pulmonary;

        model Parameters_Systemic
          parameter Real C_svl(unit = "m6.J-1") = 0.0037509e-6;
          parameter Real C_svn(unit = "m6.J-1") = 0.1125281e-6;
          parameter Real C_svc(unit = "m6.J-1") = 0.0375094e-6;
          parameter Real C_ivl(unit = "m6.J-1") = 0.0112528e-6;
          parameter Real C_ivn(unit = "m6.J-1") = 0.5626407e-6;
          parameter Real C_ivc(unit = "m6.J-1") = 0.1125281e-6;
          parameter Real R_svl(unit = "J.s.m-6") = 18.662e+6;
          parameter Real R_svn(unit = "J.s.m-6") = 3.999e+6;
          parameter Real R_svc(unit = "J.s.m-6") = 0.06665e+6;
          parameter Real R_ivl(unit = "J.s.m-6") = 5.332e+6;
          parameter Real R_ivn(unit = "J.s.m-6") = 1.1997e+6;
          parameter Real R_ivc(unit = "J.s.m-6") = 0.06665e+6;
          parameter Real I_svl(unit = "J.s2.m-6") = 0.1333e+6;
          parameter Real I_svn(unit = "J.s2.m-6") = 0.06665e+6;
          parameter Real I_svc(unit = "J.s2.m-6") = 0.06665e+6;
          parameter Real I_ivl(unit = "J.s2.m-6") = 0.1333e+6;
          parameter Real I_ivn(unit = "J.s2.m-6") = 0.06665e+6;
          parameter Real I_ivc(unit = "J.s2.m-6") = 0.06665e+6;
          parameter Real r_ascending_aorta_A(unit = "m") = 15.642e-03;
          parameter Real r_ascending_aorta_B(unit = "m") = 15.08e-03;
          parameter Real r_ascending_aorta_C(unit = "m") = 14.5305e-03;
          parameter Real r_ascending_aorta_D(unit = "m") = 13.914e-03;
          parameter Real r_aortic_arch_C2(unit = "m") = 13.3364e-03;
          parameter Real r_brachiocephalic_trunk_C4(unit = "m") = 6.41887e-03;
          parameter Real r_aortic_arch_C46(unit = "m") = 12.76710e-03;
          parameter Real r_aortic_arch_C64(unit = "m") = 12.42880e-03;
          parameter Real r_aortic_arch_C94(unit = "m") = 11.7401e-03;
          parameter Real r_thoracic_aorta_C96(unit = "m") = 10.4579e-03;
          parameter Real r_thoracic_aorta_C100(unit = "m") = 10.2897e-03;
          parameter Real r_thoracic_aorta_C104(unit = "m") = 10.0681e-03;
          parameter Real r_thoracic_aorta_C108(unit = "m") = 9.87279e-03;
          parameter Real r_thoracic_aorta_C112(unit = "m") = 8.47543e-03;
          parameter Real r_abdominal_aorta_C114(unit = "m") = 7.51533e-03;
          parameter Real r_abdominal_aorta_C136(unit = "m") = 7.42666e-03;
          parameter Real r_abdominal_aorta_C164(unit = "m") = 7.29199e-03;
          parameter Real r_abdominal_aorta_C176(unit = "m") = 7.18902e-03;
          parameter Real r_abdominal_aorta_C188(unit = "m") = 6.79934e-03;
          parameter Real r_abdominal_aorta_C192(unit = "m") = 6.12422e-03;
          parameter Real r_posterior_intercostal_T1_R98(unit = "m") = 1.4e-03;
          parameter Real r_posterior_intercostal_T1_L102(unit = "m") = 1.4e-03;
          parameter Real r_posterior_intercostal_T2_R106(unit = "m") = 1.55e-03;
          parameter Real r_posterior_intercostal_T2_L110(unit = "m") = 1.55e-03;
          parameter Real r_celiac_trunk_C116(unit = "m") = 3.29653e-03;
          parameter Real r_splenic_T2_C118(unit = "m") = 2.16682e-03;
          parameter Real r_left_gastric_T3_C120(unit = "m") = 1.50666e-03;
          parameter Real r_splenic_T2_C122(unit = "m") = 2.16682e-03;
          parameter Real r_dorsal_pancreatic_T1_C124(unit = "m") = 0.558491e-03;
          parameter Real r_splenic_T2_C126(unit = "m") = 2.16682e-03;
          parameter Real r_common_hepatic_C128(unit = "m") = 2.68614e-03;
          parameter Real r_hepatic_artery_proper_C130(unit = "m") = 1.77555e-03;
          parameter Real r_hepatic_artery_proper_left_branch_C132(unit = "m") = 1.1663e-03;
          parameter Real r_hepatic_artery_proper_right_branch_C134(unit = "m") = 1.42068e-03;
          parameter Real r_superior_mesenteric_T4_C138(unit = "m") = 3.72737e-03;
          parameter Real r_middle_colic_T8_C140(unit = "m") = 1.425e-03;
          parameter Real r_superior_mesenteric_T4_C142(unit = "m") = 3.40146e-03;
          parameter Real r_jejunal_3_T10_C144(unit = "m") = 1.58037e-03;
          parameter Real r_superior_mesenteric_T4_C146(unit = "m") = 3.06914e-03;
          parameter Real r_jejunal_6_T11_C148(unit = "m") = 1.58037e-03;
          parameter Real r_superior_mesenteric_T4_C150(unit = "m") = 2.85201e-03;
          parameter Real r_ileocolic_T9_C152(unit = "m") = 2.0e-03;
          parameter Real r_superior_mesenteric_T4_C154(unit = "m") = 2.69257e-03;
          parameter Real r_ileal_4_T12_C156(unit = "m") = 1.8015e-03;
          parameter Real r_superior_mesenteric_T4_C158(unit = "m") = 2.49385e-03;
          parameter Real r_ileal_6_T13_C160(unit = "m") = 1.80150e-03;
          parameter Real r_superior_mesenteric_T4_C162(unit = "m") = 2.18233e-03;
          parameter Real r_renal_L166(unit = "m") = 2.73574e-03;
          parameter Real r_renal_anterior_branch_L168(unit = "m") = 2.48193e-03;
          parameter Real r_inferior_segmental_T5_L170(unit = "m") = 1.92732e-03;
          parameter Real r_superior_segmental_T4_L172(unit = "m") = 1.92732e-03;
          parameter Real r_renal_posterior_branch_T3_L174(unit = "m") = 1.59319e-03;
          parameter Real r_renal_R178(unit = "m") = 2.96767e-03;
          parameter Real r_renal_anterior_branch_R180(unit = "m") = 2.48193e-03;
          parameter Real r_superior_segmental_T4_R182(unit = "m") = 1.92732e-03;
          parameter Real r_inferior_segmental_T5_R184(unit = "m") = 1.92732e-03;
          parameter Real r_renal_posterior_branch_T3_R186(unit = "m") = 1.59319e-03;
          parameter Real r_inferior_mesenteric_T5_C190(unit = "m") = 2.07748e-03;
          parameter Real r_common_iliac_R216(unit = "m") = 4.30633e-03;
          parameter Real r_internal_iliac_T1_R218(unit = "m") = 2.81829e-03;
          parameter Real r_external_iliac_R220(unit = "m") = 3.28821e-03;
          parameter Real r_femoral_R222(unit = "m") = 3.17347e-03;
          parameter Real r_profundus_T2_R224(unit = "m") = 2.14445e-03;
          parameter Real r_femoral_R226(unit = "m") = 2.89103e-03;
          parameter Real r_popliteal_R228(unit = "m") = 2.51554e-03;
          parameter Real r_anterior_tibial_T3_R230(unit = "m") = 1.1663e-03;
          parameter Real r_popliteal_R232(unit = "m") = 2.35852e-03;
          parameter Real r_tibiofibular_trunk_R234(unit = "m") = 2.34646e-03;
          parameter Real r_posterior_tibial_T4_R236(unit = "m") = 1.22936e-03;
          parameter Real r_common_iliac_L194(unit = "m") = 4.28142e-03;
          parameter Real r_internal_iliac_T1_L196(unit = "m") = 2.81829e-03;
          parameter Real r_external_iliac_L198(unit = "m") = 3.28821e-03;
          parameter Real r_femoral_L200(unit = "m") = 3.17347e-03;
          parameter Real r_profundus_T2_L202(unit = "m") = 2.14445e-03;
          parameter Real r_femoral_L204(unit = "m") = 2.89103e-03;
          parameter Real r_popliteal_L206(unit = "m") = 2.51554e-03;
          parameter Real r_anterior_tibial_T3_L208(unit = "m") = 1.1663e-03;
          parameter Real r_popliteal_L210(unit = "m") = 2.35852e-03;
          parameter Real r_tibiofibular_trunk_L212(unit = "m") = 2.34646e-03;
          parameter Real r_posterior_tibial_T4_L214(unit = "m") = 1.22936e-03;
          parameter Real r_subclavian_R28(unit = "m") = 4.52027e-03;
          parameter Real r_subclavian_R30(unit = "m") = 3.32268e-03;
          parameter Real r_axillary_R32(unit = "m") = 2.18463e-03;
          parameter Real r_brachial_R34(unit = "m") = 1.96732e-03;
          parameter Real r_ulnar_T2_R36(unit = "m") = 1.408e-03;
          parameter Real r_common_interosseous_R38(unit = "m") = 0.959006e-03;
          parameter Real r_posterior_interosseous_T3_R40(unit = "m") = 0.675992e-03;
          parameter Real r_ulnar_T2_R42(unit = "m") = 1.408e-03;
          parameter Real r_radial_T1_R44(unit = "m") = 1.378e-03;
          parameter Real r_subclavian_L66(unit = "m") = 3.99235e-03;
          parameter Real r_subclavian_L78(unit = "m") = 2.90824e-03;
          parameter Real r_axillary_L80(unit = "m") = 2.18463e-03;
          parameter Real r_brachial_L82(unit = "m") = 1.96732e-03;
          parameter Real r_ulnar_T2_L84(unit = "m") = 1.408e-03;
          parameter Real r_common_interosseous_L86(unit = "m") = 0.959006e-03;
          parameter Real r_posterior_interosseous_T3_L88(unit = "m") = 0.675992e-03;
          parameter Real r_ulnar_T2_L90(unit = "m") = 1.408e-03;
          parameter Real r_radial_T1_L92(unit = "m") = 1.378e-03;
          parameter Real r_common_carotid_R6_A(unit = "m") = 4.43053e-03;
          parameter Real r_common_carotid_R6_B(unit = "m") = 4.137e-03;
          parameter Real r_common_carotid_R6_C(unit = "m") = 3.64938e-03;
          parameter Real r_internal_carotid_R8_A(unit = "m") = 2.53763e-03;
          parameter Real r_internal_carotid_R8_B(unit = "m") = 2.04793e-03;
          parameter Real r_internal_carotid_R8_C(unit = "m") = 1.56726e-03;
          parameter Real r_external_carotid_T2_R26(unit = "m") = 2.26547e-03;
          parameter Real r_common_carotid_L48_A(unit = "m") = 4.36635e-03;
          parameter Real r_common_carotid_L48_B(unit = "m") = 4.12756e-03;
          parameter Real r_common_carotid_L48_C(unit = "m") = 3.92047e-03;
          parameter Real r_common_carotid_L48_D(unit = "m") = 3.57978e-03;
          parameter Real r_internal_carotid_L50_A(unit = "m") = 2.53763e-03;
          parameter Real r_internal_carotid_L50_B(unit = "m") = 2.04793e-03;
          parameter Real r_internal_carotid_L50_C(unit = "m") = 1.56726e-03;
          parameter Real r_external_carotid_T2_L62(unit = "m") = 2.26547e-03;
          parameter Real r_vertebral_L2(unit = "m") = 0.133527e-2;
          parameter Real r_vertebral_R272(unit = "m") = 0.133527e-2;
          parameter Real r_basilar_C4(unit = "m") = 0.172394e-2;
          parameter Real r_posterior_cerebral_precommunicating_part_L6(unit = "m") = 0.0816459e-2;
          parameter Real r_posterior_cerebral_precommunicating_part_R204(unit = "m") = 0.0816459e-2;
          parameter Real r_posterior_communicating_L8(unit = "m") = 0.0494463e-2;
          parameter Real r_posterior_communicating_R206(unit = "m") = 0.0494463e-2;
          parameter Real r_posterior_cerebral_postcommunicating_part_L12(unit = "m") = 0.0869072e-2;
          parameter Real r_posterior_cerebral_postcommunicating_part_R208(unit = "m") = 0.0869072e-2;
          parameter Real r_occipital_lateral_L14(unit = "m") = 0.064999e-2;
          parameter Real r_occipital_lateral_R210(unit = "m") = 0.0649994e-2;
          parameter Real r_lateral_occipital_anterior_temporal_branch_T46_L16(unit = "m") = 0.06e-2;
          parameter Real r_lateral_occipital_anterior_temporal_branch_T111_R212(unit = "m") = 0.06e-2;
          parameter Real r_occipital_lateral_L18(unit = "m") = 0.0649994e-2;
          parameter Real r_occipital_lateral_R214(unit = "m") = 0.0649994e-2;
          parameter Real r_lateral_occipital_intermediate_temporal_branch_T47_L20(unit = "m") = 0.028e-2;
          parameter Real r_lateral_occipital_intermediate_temporal_branch_T76_R216(unit = "m") = 0.028e-2;
          parameter Real r_occipital_lateral_L22(unit = "m") = 0.0649994e-2;
          parameter Real r_occipital_lateral_R218(unit = "m") = 0.0649994e-2;
          parameter Real r_lateral_occipital_posterior_temporal_branch_T48_L24(unit = "m") = 0.0375006e-2;
          parameter Real r_lateral_occipital_posterior_temporal_branch_T112_R220(unit = "m") = 0.0375006e-2;
          parameter Real r_medial_occipital_L26(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_R222(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_dorsal_branch_to_corpus_callosum_L28(unit = "m") = 0.0349979e-2;
          parameter Real r_medial_occipital_dorsal_branch_to_corpus_callosum_R224(unit = "m") = 0.0349979e-2;
          parameter Real r_pericallosal_parieto_occipital_branch_T60_L30(unit = "m") = 0.0421522e-2;
          parameter Real r_pericallosal_parieto_occipital_branch_T124_R226(unit = "m") = 0.0421522e-2;
          parameter Real r_pericallosal_L32(unit = "m") = 0.0800314e-2;
          parameter Real r_pericallosal_R228(unit = "m") = 0.0800314e-2;
          parameter Real r_pericallosal_precuneal_branch_T61_L34(unit = "m") = 0.0477132e-2;
          parameter Real r_pericallosal_precuneal_branch_T125_R230(unit = "m") = 0.0477132e-2;
          parameter Real r_pericallosal_L36(unit = "m") = 0.0800314e-2;
          parameter Real r_pericallosal_R232(unit = "m") = 0.0800314e-2;
          parameter Real r_anterior_cerebral_L38(unit = "m") = 0.0965045e-2;
          parameter Real r_anterior_cerebral_R234(unit = "m") = 0.0965045e-2;
          parameter Real r_distal_medial_striate_T44_L40(unit = "m") = 0.027251e-2;
          parameter Real r_distal_medial_striate_T109_R236(unit = "m") = 0.027251e-2;
          parameter Real r_anterior_cerebral_L42(unit = "m") = 0.0965045e-2;
          parameter Real r_anterior_cerebral_R238(unit = "m") = 0.0965045e-2;
          parameter Real r_anterior_communicating_C44(unit = "m") = 0.0656645e-2;
          parameter Real r_anterior_cerebral_L110(unit = "m") = 0.0965045e-2;
          parameter Real r_anterior_cerebral_R46(unit = "m") = 0.0965045e-2;
          parameter Real r_internal_carotid_L112(unit = "m") = 0.132843e-2;
          parameter Real r_internal_carotid_R48(unit = "m") = 0.132843e-2;
          parameter Real r_middle_cerebral_L114(unit = "m") = 0.104132e-2;
          parameter Real r_middle_cerebral_R52(unit = "m") = 0.104132e-2;
          parameter Real r_anterior_choroidal_T34_L116(unit = "m") = 0.0519606e-2;
          parameter Real r_anterior_choroidal_T98_R54(unit = "m") = 0.0519606e-2;
          parameter Real r_middle_cerebral_L118(unit = "m") = 0.104132e-2;
          parameter Real r_middle_cerebral_R56(unit = "m") = 0.104132e-2;
          parameter Real r_middle_cerebral_superior_terminal_branch_L120(unit = "m") = 0.0877201e-2;
          parameter Real r_middle_cerebral_superior_terminal_branch_R58(unit = "m") = 0.0877201e-2;
          parameter Real r_lateral_frontobasal_T45_L122(unit = "m") = 0.0385346e-2;
          parameter Real r_lateral_frontobasal_T110_R60(unit = "m") = 0.0385346e-2;
          parameter Real r_middle_cerebral_superior_terminal_branch_L124(unit = "m") = 0.0877201e-2;
          parameter Real r_middle_cerebral_superior_terminal_branch_R62(unit = "m") = 0.0877201e-2;
          parameter Real r_prefrontal_T65_L126(unit = "m") = 0.0480986e-2;
          parameter Real r_prefrontal_T130_R64(unit = "m") = 0.0480986e-2;
          parameter Real r_middle_cerebral_superior_terminal_branch_L128(unit = "m") = 0.0877201e-2;
          parameter Real r_middle_cerebral_superior_terminal_branch_R66(unit = "m") = 0.0877201e-2;
          parameter Real r_artery_of_precentral_sulcus_T38_L130(unit = "m") = 0.0480986e-2;
          parameter Real r_artery_of_precentral_sulcus_T103_R68(unit = "m") = 0.0480986e-2;
          parameter Real r_artery_of_central_sulcus_T36_L132(unit = "m") = 0.0480986e-2;
          parameter Real r_artery_of_central_sulcus_T101_R70(unit = "m") = 0.0480986e-2;
          parameter Real r_artery_of_precentral_sulcus_T38_L134(unit = "m") = 0.0480986e-2;
          parameter Real r_artery_of_precentral_sulcus_T103_R72(unit = "m") = 0.0480986e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_L136(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_R74(unit = "m") = 0.0816459e-2;
          parameter Real r_polar_temporal_T63_L138(unit = "m") = 0.0309894e-2;
          parameter Real r_polar_temporal_T127_R76(unit = "m") = 0.0309894e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_L140(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_R78(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_anterior_temporal_branch_T55_L142(unit = "m") = 0.0442556e-2;
          parameter Real r_middle_cerebral_anterior_temporal_branch_T119_R80(unit = "m") = 0.0442556e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_L144(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_R82(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_middle_temporal_branch_T57_L146(unit = "m") = 0.046174e-2;
          parameter Real r_middle_cerebral_middle_temporal_branch_T121_R84(unit = "m") = 0.046174e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_L148(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_R86(unit = "m") = 0.0816459e-2;
          parameter Real r_artery_of_postcentral_sulcus_T37_L150(unit = "m") = 0.0480986e-2;
          parameter Real r_artery_of_postcentral_sulcus_T102_R88(unit = "m") = 0.0480986e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_L152(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_R90(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_posterior_temporal_branch_T58_L154(unit = "m") = 0.0500255e-2;
          parameter Real r_middle_cerebral_posterior_temporal_branch_T122_R92(unit = "m") = 0.0500255e-2;
          parameter Real r_middle_cerebral_posterior_temporo_occipital_branch_T59_L156(unit = "m") = 0.0500255e-2;
          parameter Real r_middle_cerebral_posterior_temporo_occipital_branch_T123_R94(unit = "m") = 0.0500255e-2;
          parameter Real r_middle_cerebral_posterior_temporal_branch_T58_L158(unit = "m") = 0.0500255e-2;
          parameter Real r_middle_cerebral_posterior_temporal_branch_T122_R96(unit = "m") = 0.0500255e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_L160(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_R98(unit = "m") = 0.0816459e-2;
          parameter Real r_anterior_parietal_T35_L162(unit = "m") = 0.0442556e-2;
          parameter Real r_anterior_parietal_T100_R100(unit = "m") = 0.0442556e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_L164(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_inferior_terminal_branch_R102(unit = "m") = 0.0816459e-2;
          parameter Real r_middle_cerebral_branch_to_angular_gyrus_T56_L166(unit = "m") = 0.0558491e-2;
          parameter Real r_middle_cerebral_branch_to_angular_gyrus_T120_R104(unit = "m") = 0.0558491e-2;
          parameter Real r_posterior_parietal_T64_L168(unit = "m") = 0.0519606e-2;
          parameter Real r_posterior_parietal_T129_R106(unit = "m") = 0.0519606e-2;
          parameter Real r_middle_cerebral_branch_to_angular_gyrus_T56_L170(unit = "m") = 0.0558491e-2;
          parameter Real r_middle_cerebral_branch_to_angular_gyrus_T120_R108(unit = "m") = 0.0558491e-2;
          parameter Real r_callosomarginal_L172(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_R240(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_intermediomedial_frontal_branch_T41_L174(unit = "m") = 0.0465585e-2;
          parameter Real r_callosomarginal_intermediomedial_frontal_branch_T106_R242(unit = "m") = 0.0465585e-2;
          parameter Real r_callosomarginal_L176(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_R244(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_posteromedial_frontal_branch_T43_L178(unit = "m") = 0.0521563e-2;
          parameter Real r_callosomarginal_posteromedial_frontal_branch_T108_R246(unit = "m") = 0.0521563e-2;
          parameter Real r_callosomarginal_L180(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_R248(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_cingular_branch_T40_L182(unit = "m") = 0.0475194e-2;
          parameter Real r_callosomarginal_cingular_branch_T105_R250(unit = "m") = 0.0475194e-2;
          parameter Real r_callosomarginal_L184(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_R252(unit = "m") = 0.0623142e-2;
          parameter Real r_callosomarginal_paracentral_branch_T42_L186(unit = "m") = 0.0475194e-2;
          parameter Real r_callosomarginal_paracentral_branch_T107_R254(unit = "m") = 0.0475194e-2;
          parameter Real r_medial_occipital_L188(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_R256(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_occipitotemporal_branch_T52_L190(unit = "m") = 0.0349979e-2;
          parameter Real r_medial_occipital_occipitotemporal_branch_T80_R258(unit = "m") = 0.0349979e-2;
          parameter Real r_medial_occipital_L192(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_R260(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_parieto_occipital_branch_T54_L194(unit = "m") = 0.0500001e-2;
          parameter Real r_medial_occipital_parieto_occipital_branch_T118_R262(unit = "m") = 0.0500001e-2;
          parameter Real r_medial_occipital_parietal_branch_T53_L196(unit = "m") = 0.0295002e-2;
          parameter Real r_medial_occipital_parietal_branch_T81_R264(unit = "m") = 0.0295002e-2;
          parameter Real r_medial_occipital_parieto_occipital_branch_T54_L198(unit = "m") = 0.0500001e-2;
          parameter Real r_medial_occipital_parieto_occipital_branch_T118_R266(unit = "m") = 0.0500001e-2;
          parameter Real r_medial_occipital_L200(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_R268(unit = "m") = 0.0795007e-2;
          parameter Real r_medial_occipital_calcarine_branch_T51_L202(unit = "m") = 0.0349979e-2;
          parameter Real r_medial_occipital_calcarine_branch_T79_R270(unit = "m") = 0.0349979e-2;
          parameter Real l_ascending_aorta_A(unit = "m") = 15.3234e-03;
          parameter Real l_ascending_aorta_B(unit = "m") = 14.796e-03;
          parameter Real l_ascending_aorta_C(unit = "m") = 14.796e-03;
          parameter Real l_ascending_aorta_D(unit = "m") = 14.796e-03;
          parameter Real l_aortic_arch_C2(unit = "m") = 14.796e-03;
          parameter Real l_brachiocephalic_trunk_C4(unit = "m") = 47.3822e-03;
          parameter Real l_aortic_arch_C46(unit = "m") = 9.60849e-03;
          parameter Real l_aortic_arch_C64(unit = "m") = 6.97955e-03;
          parameter Real l_aortic_arch_C94(unit = "m") = 43.2111e-03;
          parameter Real l_thoracic_aorta_C96(unit = "m") = 9.89803e-03;
          parameter Real l_thoracic_aorta_C100(unit = "m") = 7.88038e-03;
          parameter Real l_thoracic_aorta_C104(unit = "m") = 15.5561e-03;
          parameter Real l_thoracic_aorta_C108(unit = "m") = 5.32705e-03;
          parameter Real l_thoracic_aorta_C112(unit = "m") = 121.566e-03;
          parameter Real l_abdominal_aorta_C114(unit = "m") = 3.24767e-03;
          parameter Real l_abdominal_aorta_C136(unit = "m") = 13.9886e-03;
          parameter Real l_abdominal_aorta_C164(unit = "m") = 4.31913e-03;
          parameter Real l_abdominal_aorta_C176(unit = "m") = 11.9773e-03;
          parameter Real l_abdominal_aorta_C188(unit = "m") = 54.0907e-03;
          parameter Real l_abdominal_aorta_C192(unit = "m") = 42.231e-03;
          parameter Real l_posterior_intercostal_T1_R98(unit = "m") = 197.232e-03;
          parameter Real l_posterior_intercostal_T1_L102(unit = "m") = 178.519e-03;
          parameter Real l_posterior_intercostal_T2_R106(unit = "m") = 201.883e-03;
          parameter Real l_posterior_intercostal_T2_L110(unit = "m") = 185.547e-03;
          parameter Real l_celiac_trunk_C116(unit = "m") = 16.9374e-03;
          parameter Real l_splenic_T2_C118(unit = "m") = 3.9576e-03;
          parameter Real l_left_gastric_T3_C120(unit = "m") = 94.8344e-03;
          parameter Real l_splenic_T2_C122(unit = "m") = 2.79812e-03;
          parameter Real l_dorsal_pancreatic_T1_C124(unit = "m") = 33.4687e-03;
          parameter Real l_splenic_T2_C126(unit = "m") = 63.2749e-03;
          parameter Real l_common_hepatic_C128(unit = "m") = 69.3076e-03;
          parameter Real l_hepatic_artery_proper_C130(unit = "m") = 16.8059e-03;
          parameter Real l_hepatic_artery_proper_left_branch_C132(unit = "m") = 164.224e-03;
          parameter Real l_hepatic_artery_proper_right_branch_C134(unit = "m") = 80.0632e-03;
          parameter Real l_superior_mesenteric_T4_C138(unit = "m") = 49.5492e-03;
          parameter Real l_middle_colic_T8_C140(unit = "m") = 116.298e-03;
          parameter Real l_superior_mesenteric_T4_C142(unit = "m") = 35.1664e-03;
          parameter Real l_jejunal_3_T10_C144(unit = "m") = 47.6944e-03;
          parameter Real l_superior_mesenteric_T4_C146(unit = "m") = 32.2488e-03;
          parameter Real l_jejunal_6_T11_C148(unit = "m") = 63.8535e-03;
          parameter Real l_superior_mesenteric_T4_C150(unit = "m") = 16.7458e-03;
          parameter Real l_ileocolic_T9_C152(unit = "m") = 46.8332e-03;
          parameter Real l_superior_mesenteric_T4_C154(unit = "m") = 23.0546e-03;
          parameter Real l_ileal_4_T12_C156(unit = "m") = 45.6984e-03;
          parameter Real l_superior_mesenteric_T4_C158(unit = "m") = 20.5397e-03;
          parameter Real l_ileal_6_T13_C160(unit = "m") = 29.1584e-03;
          parameter Real l_superior_mesenteric_T4_C162(unit = "m") = 39.4879e-03;
          parameter Real l_renal_L166(unit = "m") = 22.0037e-03;
          parameter Real l_renal_anterior_branch_L168(unit = "m") = 10.8789e-03;
          parameter Real l_inferior_segmental_T5_L170(unit = "m") = 40.8761e-03;
          parameter Real l_superior_segmental_T4_L172(unit = "m") = 29.7265e-03;
          parameter Real l_renal_posterior_branch_T3_L174(unit = "m") = 22.3608e-03;
          parameter Real l_renal_R178(unit = "m") = 37.7403e-03;
          parameter Real l_renal_anterior_branch_R180(unit = "m") = 10.8792e-03;
          parameter Real l_superior_segmental_T4_R182(unit = "m") = 29.7263e-03;
          parameter Real l_inferior_segmental_T5_R184(unit = "m") = 40.8756e-03;
          parameter Real l_renal_posterior_branch_T3_R186(unit = "m") = 22.36e-03;
          parameter Real l_inferior_mesenteric_T5_C190(unit = "m") = 90.3282e-03;
          parameter Real l_common_iliac_R216(unit = "m") = 76.4393e-03;
          parameter Real l_internal_iliac_T1_R218(unit = "m") = 72.5302e-03;
          parameter Real l_external_iliac_R220(unit = "m") = 102.358e-03;
          parameter Real l_femoral_R222(unit = "m") = 31.5982e-03;
          parameter Real l_profundus_T2_R224(unit = "m") = 238.438e-03;
          parameter Real l_femoral_R226(unit = "m") = 319.297e-03;
          parameter Real l_popliteal_R228(unit = "m") = 132.06e-03;
          parameter Real l_anterior_tibial_T3_R230(unit = "m") = 386.388e-03;
          parameter Real l_popliteal_R232(unit = "m") = 8.80051e-03;
          parameter Real l_tibiofibular_trunk_R234(unit = "m") = 36.1667e-03;
          parameter Real l_posterior_tibial_T4_R236(unit = "m") = 382.987e-03;
          parameter Real l_common_iliac_L194(unit = "m") = 74.0524e-03;
          parameter Real l_internal_iliac_T1_L196(unit = "m") = 72.5301e-03;
          parameter Real l_external_iliac_L198(unit = "m") = 102.358e-03;
          parameter Real l_femoral_L200(unit = "m") = 31.5982e-03;
          parameter Real l_profundus_T2_L202(unit = "m") = 238.438e-03;
          parameter Real l_femoral_L204(unit = "m") = 319.297e-03;
          parameter Real l_popliteal_L206(unit = "m") = 132.059e-03;
          parameter Real l_anterior_tibial_T3_L208(unit = "m") = 386.389e-03;
          parameter Real l_popliteal_L210(unit = "m") = 8.80046e-03;
          parameter Real l_tibiofibular_trunk_L212(unit = "m") = 36.1676e-03;
          parameter Real l_posterior_tibial_T4_L214(unit = "m") = 382.987e-03;
          parameter Real l_subclavian_R28(unit = "m") = 15.7469e-03;
          parameter Real l_subclavian_R30(unit = "m") = 41.1419e-03;
          parameter Real l_axillary_R32(unit = "m") = 120.021e-03;
          parameter Real l_brachial_R34(unit = "m") = 223.119e-03;
          parameter Real l_ulnar_T2_R36(unit = "m") = 29.7599e-03;
          parameter Real l_common_interosseous_R38(unit = "m") = 16.2682e-03;
          parameter Real l_posterior_interosseous_T3_R40(unit = "m") = 231.694e-03;
          parameter Real l_ulnar_T2_R42(unit = "m") = 239.276e-03;
          parameter Real l_radial_T1_R44(unit = "m") = 302.156e-03;
          parameter Real l_subclavian_L66(unit = "m") = 49.4669e-03;
          parameter Real l_subclavian_L78(unit = "m") = 41.1396e-03;
          parameter Real l_axillary_L80(unit = "m") = 120.021e-03;
          parameter Real l_brachial_L82(unit = "m") = 223.119e-03;
          parameter Real l_ulnar_T2_L84(unit = "m") = 29.7594e-03;
          parameter Real l_common_interosseous_L86(unit = "m") = 16.2681e-03;
          parameter Real l_posterior_interosseous_T3_L88(unit = "m") = 231.695e-03;
          parameter Real l_ulnar_T2_L90(unit = "m") = 239.277e-03;
          parameter Real l_radial_T1_L92(unit = "m") = 302.155e-03;
          parameter Real l_common_carotid_R6_A(unit = "m") = 27.0844e-03;
          parameter Real l_common_carotid_R6_B(unit = "m") = 27.0844e-03;
          parameter Real l_common_carotid_R6_C(unit = "m") = 27.0844e-03;
          parameter Real l_internal_carotid_R8_A(unit = "m") = 45.036e-03;
          parameter Real l_internal_carotid_R8_B(unit = "m") = 45.036e-03;
          parameter Real l_internal_carotid_R8_C(unit = "m") = 45.036e-03;
          parameter Real l_external_carotid_T2_R26(unit = "m") = 61.0125e-03;
          parameter Real l_common_carotid_L48_A(unit = "m") = 30.339e-03;
          parameter Real l_common_carotid_L48_B(unit = "m") = 30.339e-03;
          parameter Real l_common_carotid_L48_C(unit = "m") = 30.339e-03;
          parameter Real l_common_carotid_L48_D(unit = "m") = 30.339e-03;
          parameter Real l_internal_carotid_L50_A(unit = "m") = 45.036e-03;
          parameter Real l_internal_carotid_L50_B(unit = "m") = 45.036e-03;
          parameter Real l_internal_carotid_L50_C(unit = "m") = 45.036e-03;
          parameter Real l_external_carotid_T2_L62(unit = "m") = 61.0127e-03;
          parameter Real l_vertebral_L2(unit = "m") = 20.9765e-2;
          parameter Real l_vertebral_R272(unit = "m") = 21.0146e-2;
          parameter Real l_basilar_C4(unit = "m") = 2.26443e-2;
          parameter Real l_posterior_cerebral_precommunicating_part_L6(unit = "m") = 0.697753e-2;
          parameter Real l_posterior_cerebral_precommunicating_part_R204(unit = "m") = 0.77827e-2;
          parameter Real l_posterior_communicating_L8(unit = "m") = 1.65974e-2;
          parameter Real l_posterior_communicating_R206(unit = "m") = 1.6597e-2;
          parameter Real l_posterior_cerebral_postcommunicating_part_L12(unit = "m") = 1.28091e-2;
          parameter Real l_posterior_cerebral_postcommunicating_part_R208(unit = "m") = 1.28083e-2;
          parameter Real l_occipital_lateral_L14(unit = "m") = 0.472985e-2;
          parameter Real l_occipital_lateral_R210(unit = "m") = 0.472973e-2;
          parameter Real l_lateral_occipital_anterior_temporal_branch_T111_R212(unit = "m") = 2.61986e-2;
          parameter Real l_lateral_occipital_anterior_temporal_branch_T46_L16(unit = "m") = 2.61979e-2;
          parameter Real l_occipital_lateral_L18(unit = "m") = 1.11052e-2;
          parameter Real l_occipital_lateral_R214(unit = "m") = 1.11056e-2;
          parameter Real l_lateral_occipital_intermediate_temporal_branch_T47_L20(unit = "m") = 2.1112e-2;
          parameter Real l_lateral_occipital_intermediate_temporal_branch_T76_R216(unit = "m") = 2.1112e-2;
          parameter Real l_occipital_lateral_L22(unit = "m") = 1.2455e-2;
          parameter Real l_occipital_lateral_R218(unit = "m") = 1.2455e-2;
          parameter Real l_lateral_occipital_posterior_temporal_branch_T48_L24(unit = "m") = 2.05518e-2;
          parameter Real l_lateral_occipital_posterior_temporal_branch_T112_R220(unit = "m") = 2.05518e-2;
          parameter Real l_medial_occipital_L26(unit = "m") = 3.21752e-2;
          parameter Real l_medial_occipital_R222(unit = "m") = 3.21754e-2;
          parameter Real l_medial_occipital_dorsal_branch_to_corpus_callosum_L28(unit = "m") = 2.51851e-2;
          parameter Real l_medial_occipital_dorsal_branch_to_corpus_callosum_R224(unit = "m") = 2.51847e-2;
          parameter Real l_pericallosal_parieto_occipital_branch_T60_L30(unit = "m") = 4.1888e-2;
          parameter Real l_pericallosal_parieto_occipital_branch_T124_R226(unit = "m") = 4.18881e-2;
          parameter Real l_pericallosal_L32(unit = "m") = 1.5668e-2;
          parameter Real l_pericallosal_R228(unit = "m") = 1.5668e-2;
          parameter Real l_pericallosal_precuneal_branch_T61_L34(unit = "m") = 5.42988e-2;
          parameter Real l_pericallosal_precuneal_branch_T125_R230(unit = "m") = 5.42989e-2;
          parameter Real l_pericallosal_L36(unit = "m") = 5.57509e-2;
          parameter Real l_pericallosal_R232(unit = "m") = 5.57509e-2;
          parameter Real l_anterior_cerebral_L38(unit = "m") = 2.809e-2;
          parameter Real l_anterior_cerebral_R234(unit = "m") = 2.80898e-2;
          parameter Real l_distal_medial_striate_T44_L40(unit = "m") = 1.24626e-2;
          parameter Real l_distal_medial_striate_T109_R236(unit = "m") = 1.24626e-2;
          parameter Real l_anterior_cerebral_L42(unit = "m") = 0.114668e-2;
          parameter Real l_anterior_cerebral_R238(unit = "m") = 0.114733e-2;
          parameter Real l_anterior_communicating_C44(unit = "m") = 0.565531e-2;
          parameter Real l_anterior_cerebral_L110(unit = "m") = 1.05951e-2;
          parameter Real l_anterior_cerebral_R46(unit = "m") = 1.05945e-2;
          parameter Real l_internal_carotid_L112(unit = "m") = 0.196966e-2;
          parameter Real l_internal_carotid_R48(unit = "m") = 0.196955e-2;
          parameter Real l_middle_cerebral_L114(unit = "m") = 0.182233e-2;
          parameter Real l_middle_cerebral_R52(unit = "m") = 0.182229e-2;
          parameter Real l_anterior_choroidal_T34_L116(unit = "m") = 2.5869e-2;
          parameter Real l_anterior_choroidal_T98_R54(unit = "m") = 2.5869e-2;
          parameter Real l_middle_cerebral_L118(unit = "m") = 2.83007e-2;
          parameter Real l_middle_cerebral_R56(unit = "m") = 2.8301e-2;
          parameter Real l_middle_cerebral_superior_terminal_branch_L120(unit = "m") = 0.251165e-2;
          parameter Real l_middle_cerebral_superior_terminal_branch_R58(unit = "m") = 0.251165e-2;
          parameter Real l_lateral_frontobasal_T45_L122(unit = "m") = 4.32804e-2;
          parameter Real l_lateral_frontobasal_T110_R60(unit = "m") = 4.32805e-2;
          parameter Real l_middle_cerebral_superior_terminal_branch_L124(unit = "m") = 1.33091e-2;
          parameter Real l_middle_cerebral_superior_terminal_branch_R62(unit = "m") = 1.3309e-2;
          parameter Real l_prefrontal_T65_L126(unit = "m") = 10.0723e-2;
          parameter Real l_prefrontal_T130_R64(unit = "m") = 10.0723e-2;
          parameter Real l_middle_cerebral_superior_terminal_branch_L128(unit = "m") = 1.2255e-2;
          parameter Real l_middle_cerebral_superior_terminal_branch_R66(unit = "m") = 1.2255e-2;
          parameter Real l_artery_of_precentral_sulcus_T38_L130(unit = "m") = 0.935653e-2;
          parameter Real l_artery_of_precentral_sulcus_T103_R68(unit = "m") = 0.935654e-2;
          parameter Real l_artery_of_central_sulcus_T36_L132(unit = "m") = 10.6498e-2;
          parameter Real l_artery_of_central_sulcus_T101_R70(unit = "m") = 10.6498e-2;
          parameter Real l_artery_of_precentral_sulcus_T38_L134(unit = "m") = 8.81918e-2;
          parameter Real l_artery_of_precentral_sulcus_T103_R72(unit = "m") = 8.81919e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_L136(unit = "m") = 0.142512e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_R74(unit = "m") = 0.142516e-2;
          parameter Real l_polar_temporal_T63_L138(unit = "m") = 2.17951e-2;
          parameter Real l_polar_temporal_T127_R76(unit = "m") = 2.17953e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_L140(unit = "m") = 0.0932844e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_R78(unit = "m") = 0.0932844e-2;
          parameter Real l_middle_cerebral_anterior_temporal_branch_T55_L142(unit = "m") = 2.67934e-2;
          parameter Real l_middle_cerebral_anterior_temporal_branch_T119_R80(unit = "m") = 2.67942e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_L144(unit = "m") = 2.83599e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_R82(unit = "m") = 2.83601e-2;
          parameter Real l_middle_cerebral_middle_temporal_branch_T57_L146(unit = "m") = 3.00383e-2;
          parameter Real l_middle_cerebral_middle_temporal_branch_T121_R84(unit = "m") = 3.00379e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_L148(unit = "m") = 1.82043e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_R86(unit = "m") = 1.82041e-2;
          parameter Real l_artery_of_postcentral_sulcus_T37_L150(unit = "m") = 7.36477e-2;
          parameter Real l_artery_of_postcentral_sulcus_T102_R88(unit = "m") = 7.36471e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_L152(unit = "m") = 0.708809e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_R90(unit = "m") = 0.708809e-2;
          parameter Real l_middle_cerebral_posterior_temporal_branch_T58_L154(unit = "m") = 0.692097e-2;
          parameter Real l_middle_cerebral_posterior_temporal_branch_T122_R92(unit = "m") = 0.692096e-2;
          parameter Real l_middle_cerebral_posterior_temporo_occipital_branch_T59_L156(unit = "m") = 5.73679e-2;
          parameter Real l_middle_cerebral_posterior_temporo_occipital_branch_T123_R94(unit = "m") = 5.7368e-2;
          parameter Real l_middle_cerebral_posterior_temporal_branch_T58_L158(unit = "m") = 2.26651e-2;
          parameter Real l_middle_cerebral_posterior_temporal_branch_T122_R96(unit = "m") = 2.26651e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_L160(unit = "m") = 0.549683e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_R98(unit = "m") = 0.549683e-2;
          parameter Real l_anterior_parietal_T35_L162(unit = "m") = 6.219e-2;
          parameter Real l_anterior_parietal_T100_R100(unit = "m") = 6.21893e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_L164(unit = "m") = 1.05444e-2;
          parameter Real l_middle_cerebral_inferior_terminal_branch_R102(unit = "m") = 1.05443e-2;
          parameter Real l_middle_cerebral_branch_to_angular_gyrus_T56_L166(unit = "m") = 1.01852e-2;
          parameter Real l_middle_cerebral_branch_to_angular_gyrus_T120_R104(unit = "m") = 1.01852e-2;
          parameter Real l_posterior_parietal_T64_L168(unit = "m") = 4.2951e-2;
          parameter Real l_posterior_parietal_T129_R106(unit = "m") = 4.2951e-2;
          parameter Real l_middle_cerebral_branch_to_angular_gyrus_T56_L170(unit = "m") = 4.19474e-2;
          parameter Real l_middle_cerebral_branch_to_angular_gyrus_T120_R108(unit = "m") = 4.19475e-2;
          parameter Real l_callosomarginal_L172(unit = "m") = 1.13312e-2;
          parameter Real l_callosomarginal_R240(unit = "m") = 1.13311e-2;
          parameter Real l_callosomarginal_intermediomedial_frontal_branch_T41_L174(unit = "m") = 3.43749e-2;
          parameter Real l_callosomarginal_intermediomedial_frontal_branch_T106_R242(unit = "m") = 3.43749e-2;
          parameter Real l_callosomarginal_L176(unit = "m") = 1.72952e-2;
          parameter Real l_callosomarginal_R244(unit = "m") = 1.72952e-2;
          parameter Real l_callosomarginal_posteromedial_frontal_branch_T43_L178(unit = "m") = 4.19271e-2;
          parameter Real l_callosomarginal_posteromedial_frontal_branch_T108_R246(unit = "m") = 4.1927e-2;
          parameter Real l_callosomarginal_L180(unit = "m") = 1.39198e-2;
          parameter Real l_callosomarginal_R248(unit = "m") = 1.39198e-2;
          parameter Real l_callosomarginal_cingular_branch_T40_L182(unit = "m") = 1.75028e-2;
          parameter Real l_callosomarginal_cingular_branch_T105_R250(unit = "m") = 1.75028e-2;
          parameter Real l_callosomarginal_L184(unit = "m") = 0.496294e-2;
          parameter Real l_callosomarginal_R252(unit = "m") = 0.496295e-2;
          parameter Real l_callosomarginal_paracentral_branch_T42_L186(unit = "m") = 4.22964e-2;
          parameter Real l_callosomarginal_paracentral_branch_T107_R254(unit = "m") = 4.22964e-2;
          parameter Real l_medial_occipital_L188(unit = "m") = 0.490507e-2;
          parameter Real l_medial_occipital_R256(unit = "m") = 0.490474e-2;
          parameter Real l_medial_occipital_occipitotemporal_branch_T52_L190(unit = "m") = 2.98359e-2;
          parameter Real l_medial_occipital_occipitotemporal_branch_T80_R258(unit = "m") = 2.98366e-2;
          parameter Real l_medial_occipital_L192(unit = "m") = 0.892062e-2;
          parameter Real l_medial_occipital_R260(unit = "m") = 0.892067e-2;
          parameter Real l_medial_occipital_parieto_occipital_branch_T54_L194(unit = "m") = 1.88987e-2;
          parameter Real l_medial_occipital_parieto_occipital_branch_T118_R262(unit = "m") = 1.88989e-2;
          parameter Real l_medial_occipital_parietal_branch_T53_L196(unit = "m") = 2.80168e-2;
          parameter Real l_medial_occipital_parietal_branch_T81_R264(unit = "m") = 2.80169e-2;
          parameter Real l_medial_occipital_parieto_occipital_branch_T54_L198(unit = "m") = 2.33801e-2;
          parameter Real l_medial_occipital_parieto_occipital_branch_T118_R266(unit = "m") = 2.33801e-2;
          parameter Real l_medial_occipital_L200(unit = "m") = 1.26049e-2;
          parameter Real l_medial_occipital_R268(unit = "m") = 1.26049e-2;
          parameter Real l_medial_occipital_calcarine_branch_T51_L202(unit = "m") = 1.47666e-2;
          parameter Real l_medial_occipital_calcarine_branch_T79_R270(unit = "m") = 1.47667e-2;
          parameter Real C_ascending_aorta_A(unit = "m6.J-1") = 1.631e-10;
          parameter Real C_ascending_aorta_B(unit = "m6.J-1") = 1.631e-10;
          parameter Real C_ascending_aorta_C(unit = "m6.J-1") = 1.631e-10;
          parameter Real C_ascending_aorta_D(unit = "m6.J-1") = 1.631e-10;
          parameter Real C_aortic_arch_C2(unit = "m6.J-1") = 1.631e-10;
          parameter Real C_brachiocephalic_trunk_C4(unit = "m6.J-1") = 0.86529e-10;
          parameter Real C_aortic_arch_C46(unit = "m6.J-1") = 0.799e-10;
          parameter Real C_aortic_arch_C64(unit = "m6.J-1") = 0.548e-10;
          parameter Real C_aortic_arch_C94(unit = "m6.J-1") = 2.83e-10;
          parameter Real C_thoracic_aorta_C96(unit = "m6.J-1") = 0.534e-10;
          parameter Real C_thoracic_aorta_C100(unit = "m6.J-1") = 0.411e-10;
          parameter Real C_thoracic_aorta_C104(unit = "m6.J-1") = 0.774e-10;
          parameter Real C_thoracic_aorta_C108(unit = "m6.J-1") = 0.254e-10;
          parameter Real C_thoracic_aorta_C112(unit = "m6.J-1") = 4.39e-10;
          parameter Real C_abdominal_aorta_C114(unit = "m6.J-1") = 0.084e-10;
          parameter Real C_abdominal_aorta_C136(unit = "m6.J-1") = 0.35e-10;
          parameter Real C_abdominal_aorta_C164(unit = "m6.J-1") = 0.104e-10;
          parameter Real C_abdominal_aorta_C176(unit = "m6.J-1") = 0.281e-10;
          parameter Real C_abdominal_aorta_C188(unit = "m6.J-1") = 1.11e-10;
          parameter Real C_abdominal_aorta_C192(unit = "m6.J-1") = 0.6966e-10;
          parameter Real C_posterior_intercostal_T1_R98(unit = "m6.J-1") = 0.0848e-10;
          parameter Real C_posterior_intercostal_T1_L102(unit = "m6.J-1") = 0.0767e-10;
          parameter Real C_posterior_intercostal_T2_R106(unit = "m6.J-1") = 0.11e-10;
          parameter Real C_posterior_intercostal_T2_L110(unit = "m6.J-1") = 0.102e-10;
          parameter Real C_celiac_trunk_C116(unit = "m6.J-1") = 0.0593e-10;
          parameter Real C_splenic_T2_C118(unit = "m6.J-1") = 0.00491e-10;
          parameter Real C_left_gastric_T3_C120(unit = "m6.J-1") = 0.04857e-10;
          parameter Real C_splenic_T2_C122(unit = "m6.J-1") = 0.00347e-10;
          parameter Real C_dorsal_pancreatic_T1_C124(unit = "m6.J-1") = 0.00179e-10;
          parameter Real C_splenic_T2_C126(unit = "m6.J-1") = 0.07847e-10;
          parameter Real C_common_hepatic_C128(unit = "m6.J-1") = 0.1468e-10;
          parameter Real C_hepatic_artery_proper_C130(unit = "m6.J-1") = 0.01254e-10;
          parameter Real C_hepatic_artery_proper_left_branch_C132(unit = "m6.J-1") = 0.046e-10;
          parameter Real C_hepatic_artery_proper_right_branch_C134(unit = "m6.J-1") = 0.0356e-10;
          parameter Real C_superior_mesenteric_T4_C138(unit = "m6.J-1") = 0.2373e-10;
          parameter Real C_middle_colic_T8_C140(unit = "m6.J-1") = 0.0521e-10;
          parameter Real C_superior_mesenteric_T4_C142(unit = "m6.J-1") = 0.13e-10;
          parameter Real C_jejunal_3_T10_C144(unit = "m6.J-1") = 0.0274e-10;
          parameter Real C_superior_mesenteric_T4_C146(unit = "m6.J-1") = 0.095e-10;
          parameter Real C_jejunal_6_T11_C148(unit = "m6.J-1") = 0.03667e-10;
          parameter Real C_superior_mesenteric_T4_C150(unit = "m6.J-1") = 0.0412e-10;
          parameter Real C_ileocolic_T9_C152(unit = "m6.J-1") = 0.04767e-10;
          parameter Real C_superior_mesenteric_T4_C154(unit = "m6.J-1") = 0.0486e-10;
          parameter Real C_ileal_4_T12_C156(unit = "m6.J-1") = 0.036e-10;
          parameter Real C_superior_mesenteric_T4_C158(unit = "m6.J-1") = 0.0361e-10;
          parameter Real C_ileal_6_T13_C160(unit = "m6.J-1") = 0.023e-10;
          parameter Real C_superior_mesenteric_T4_C162(unit = "m6.J-1") = 0.0531e-10;
          parameter Real C_renal_L166(unit = "m6.J-1") = 0.049e-10;
          parameter Real C_renal_anterior_branch_L168(unit = "m6.J-1") = 0.0189e-10;
          parameter Real C_inferior_segmental_T5_L170(unit = "m6.J-1") = 0.038e-10;
          parameter Real C_superior_segmental_T4_L172(unit = "m6.J-1") = 0.0276e-10;
          parameter Real C_renal_posterior_branch_T3_L174(unit = "m6.J-1") = 0.0131e-10;
          parameter Real C_renal_R178(unit = "m6.J-1") = 0.1e-10;
          parameter Real C_renal_anterior_branch_R180(unit = "m6.J-1") = 0.0189e-10;
          parameter Real C_superior_segmental_T4_R182(unit = "m6.J-1") = 0.02763e-10;
          parameter Real C_inferior_segmental_T5_R184(unit = "m6.J-1") = 0.038e-10;
          parameter Real C_renal_posterior_branch_T3_R186(unit = "m6.J-1") = 0.0131e-10;
          parameter Real C_inferior_mesenteric_T5_C190(unit = "m6.J-1") = 0.1e-10;
          parameter Real C_common_iliac_R216(unit = "m6.J-1") = 0.525e-10;
          parameter Real C_internal_iliac_T1_R218(unit = "m6.J-1") = 0.1733e-10;
          parameter Real C_external_iliac_R220(unit = "m6.J-1") = 0.359e-10;
          parameter Real C_femoral_R222(unit = "m6.J-1") = 0.101e-10;
          parameter Real C_profundus_T2_R224(unit = "m6.J-1") = 0.288e-10;
          parameter Real C_femoral_R226(unit = "m6.J-1") = 0.8338e-10;
          parameter Real C_popliteal_R228(unit = "m6.J-1") = 0.241e-10;
          parameter Real C_anterior_tibial_T3_R230(unit = "m6.J-1") = 0.108e-10;
          parameter Real C_popliteal_R232(unit = "m6.J-1") = 0.01345e-10;
          parameter Real C_tibiofibular_trunk_R234(unit = "m6.J-1") = 0.0546e-10;
          parameter Real C_posterior_tibial_T4_R236(unit = "m6.J-1") = 0.121e-10;
          parameter Real C_common_iliac_L194(unit = "m6.J-1") = 0.509e-10;
          parameter Real C_internal_iliac_T1_L196(unit = "m6.J-1") = 0.1733e-10;
          parameter Real C_external_iliac_L198(unit = "m6.J-1") = 0.359e-10;
          parameter Real C_femoral_L200(unit = "m6.J-1") = 0.101e-10;
          parameter Real C_profundus_T2_L202(unit = "m6.J-1") = 0.288e-10;
          parameter Real C_femoral_L204(unit = "m6.J-1") = 0.8338e-10;
          parameter Real C_popliteal_L206(unit = "m6.J-1") = 0.241e-10;
          parameter Real C_anterior_tibial_T3_L208(unit = "m6.J-1") = 0.108e-10;
          parameter Real C_popliteal_L210(unit = "m6.J-1") = 0.01345e-10;
          parameter Real C_tibiofibular_trunk_L212(unit = "m6.J-1") = 0.05466e-10;
          parameter Real C_posterior_tibial_T4_L214(unit = "m6.J-1") = 0.121e-10;
          parameter Real C_subclavian_R28(unit = "m6.J-1") = 0.1242e-10;
          parameter Real C_subclavian_R30(unit = "m6.J-1") = 0.147e-10;
          parameter Real C_axillary_R32(unit = "m6.J-1") = 0.153e-10;
          parameter Real C_brachial_R34(unit = "m6.J-1") = 0.215e-10;
          parameter Real C_ulnar_T2_R36(unit = "m6.J-1") = 0.013e-10;
          parameter Real C_common_interosseous_R38(unit = "m6.J-1") = 0.0029e-10;
          parameter Real C_posterior_interosseous_T3_R40(unit = "m6.J-1") = 0.019e-10;
          parameter Real C_ulnar_T2_R42(unit = "m6.J-1") = 0.1043e-10;
          parameter Real C_radial_T1_R44(unit = "m6.J-1") = 0.125e-10;
          parameter Real C_subclavian_L66(unit = "m6.J-1") = 0.324e-10;
          parameter Real C_subclavian_L78(unit = "m6.J-1") = 0.1075e-10;
          parameter Real C_axillary_L80(unit = "m6.J-1") = 0.153e-10;
          parameter Real C_brachial_L82(unit = "m6.J-1") = 0.215e-10;
          parameter Real C_ulnar_T2_L84(unit = "m6.J-1") = 0.013e-10;
          parameter Real C_common_interosseous_L86(unit = "m6.J-1") = 0.0029e-10;
          parameter Real C_posterior_interosseous_T3_L88(unit = "m6.J-1") = 0.019e-10;
          parameter Real C_ulnar_T2_L90(unit = "m6.J-1") = 0.1043e-10;
          parameter Real C_radial_T1_L92(unit = "m6.J-1") = 0.125e-10;
          parameter Real C_common_carotid_R6_A(unit = "m6.J-1") = 0.148e-10;
          parameter Real C_common_carotid_R6_B(unit = "m6.J-1") = 0.148e-10;
          parameter Real C_common_carotid_R6_C(unit = "m6.J-1") = 0.148e-10;
          parameter Real C_internal_carotid_R8_A(unit = "m6.J-1") = 0.0525e-10;
          parameter Real C_internal_carotid_R8_B(unit = "m6.J-1") = 0.0525e-10;
          parameter Real C_internal_carotid_R8_C(unit = "m6.J-1") = 0.0525e-10;
          parameter Real C_external_carotid_T2_R26(unit = "m6.J-1") = 0.0845e-10;
          parameter Real C_common_carotid_L48_A(unit = "m6.J-1") = 0.1661e-10;
          parameter Real C_common_carotid_L48_B(unit = "m6.J-1") = 0.1661e-10;
          parameter Real C_common_carotid_L48_C(unit = "m6.J-1") = 0.1661e-10;
          parameter Real C_common_carotid_L48_D(unit = "m6.J-1") = 0.1661e-10;
          parameter Real C_internal_carotid_L50_A(unit = "m6.J-1") = 0.0525e-10;
          parameter Real C_internal_carotid_L50_B(unit = "m6.J-1") = 0.0525e-10;
          parameter Real C_internal_carotid_L50_C(unit = "m6.J-1") = 0.0525e-10;
          parameter Real C_external_carotid_T2_L62(unit = "m6.J-1") = 0.0845e-10;
          parameter Real C_vertebral_L2(unit = "m6.J-1") = 0.0811e-10;
          parameter Real C_vertebral_R272(unit = "m6.J-1") = 0.0811e-10;
          parameter Real C_basilar_C4(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_cerebral_precommunicating_part_L6(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_cerebral_precommunicating_part_R204(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_communicating_L8(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_communicating_R206(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_cerebral_postcommunicating_part_L12(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_cerebral_postcommunicating_part_R208(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_occipital_lateral_L14(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_occipital_lateral_R210(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_occipital_anterior_temporal_branch_T111_R212(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_occipital_anterior_temporal_branch_T46_L16(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_occipital_lateral_L18(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_occipital_lateral_R214(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_occipital_intermediate_temporal_branch_T47_L20(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_occipital_intermediate_temporal_branch_T76_R216(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_occipital_lateral_L22(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_occipital_lateral_R218(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_occipital_posterior_temporal_branch_T48_L24(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_occipital_posterior_temporal_branch_T112_R220(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_L26(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_R222(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_dorsal_branch_to_corpus_callosum_L28(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_dorsal_branch_to_corpus_callosum_R224(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_parieto_occipital_branch_T60_L30(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_parieto_occipital_branch_T124_R226(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_L32(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_R228(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_precuneal_branch_T61_L34(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_precuneal_branch_T125_R230(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_L36(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_pericallosal_R232(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_cerebral_L38(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_cerebral_R234(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_distal_medial_striate_T44_L40(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_distal_medial_striate_T109_R236(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_cerebral_L42(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_cerebral_R238(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_communicating_C44(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_cerebral_L110(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_cerebral_R46(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_internal_carotid_L112(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_internal_carotid_R48(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_L114(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_R52(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_choroidal_T34_L116(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_choroidal_T98_R54(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_L118(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_R56(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_superior_terminal_branch_L120(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_superior_terminal_branch_R58(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_frontobasal_T45_L122(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_lateral_frontobasal_T110_R60(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_superior_terminal_branch_L124(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_superior_terminal_branch_R62(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_prefrontal_T65_L126(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_prefrontal_T130_R64(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_superior_terminal_branch_L128(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_superior_terminal_branch_R66(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_precentral_sulcus_T38_L130(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_precentral_sulcus_T103_R68(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_central_sulcus_T36_L132(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_central_sulcus_T101_R70(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_precentral_sulcus_T38_L134(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_precentral_sulcus_T103_R72(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_L136(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_R74(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_polar_temporal_T63_L138(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_polar_temporal_T127_R76(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_L140(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_R78(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_anterior_temporal_branch_T55_L142(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_anterior_temporal_branch_T119_R80(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_L144(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_R82(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_middle_temporal_branch_T57_L146(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_middle_temporal_branch_T121_R84(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_L148(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_R86(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_postcentral_sulcus_T37_L150(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_artery_of_postcentral_sulcus_T102_R88(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_L152(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_R90(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_posterior_temporal_branch_T58_L154(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_posterior_temporal_branch_T122_R92(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_posterior_temporo_occipital_branch_T59_L156(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_posterior_temporo_occipital_branch_T123_R94(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_posterior_temporal_branch_T58_L158(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_posterior_temporal_branch_T122_R96(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_L160(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_R98(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_parietal_T35_L162(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_anterior_parietal_T100_R100(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_L164(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_inferior_terminal_branch_R102(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_branch_to_angular_gyrus_T56_L166(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_branch_to_angular_gyrus_T120_R104(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_parietal_T64_L168(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_posterior_parietal_T129_R106(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_branch_to_angular_gyrus_T56_L170(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_middle_cerebral_branch_to_angular_gyrus_T120_R108(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_L172(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_R240(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_intermediomedial_frontal_branch_T41_L174(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_intermediomedial_frontal_branch_T106_R242(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_L176(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_R244(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_posteromedial_frontal_branch_T43_L178(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_posteromedial_frontal_branch_T108_R246(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_L180(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_R248(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_cingular_branch_T40_L182(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_cingular_branch_T105_R250(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_L184(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_R252(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_paracentral_branch_T42_L186(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_callosomarginal_paracentral_branch_T107_R254(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_L188(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_R256(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_occipitotemporal_branch_T52_L190(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_occipitotemporal_branch_T80_R258(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_L192(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_R260(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_parieto_occipital_branch_T54_L194(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_parieto_occipital_branch_T118_R262(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_parietal_branch_T53_L196(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_parietal_branch_T81_R264(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_parieto_occipital_branch_T54_L198(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_parieto_occipital_branch_T118_R266(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_L200(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_R268(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_calcarine_branch_T51_L202(unit = "m6.J-1") = 1.6e-10;
          parameter Real C_medial_occipital_calcarine_branch_T79_R270(unit = "m6.J-1") = 1.6e-10;
          parameter Real R_T_posterior_intercostal_T1_R98(unit = "J.s.m-6") = 22.5809e+9;
          parameter Real R_T_posterior_intercostal_T1_L102(unit = "J.s.m-6") = 23.1661e+9;
          parameter Real R_T_posterior_intercostal_T2_R106(unit = "J.s.m-6") = 21.0678e+9;
          parameter Real R_T_posterior_intercostal_T2_L110(unit = "J.s.m-6") = 21.2483e+9;
          parameter Real R_T_left_gastric_T3_C120(unit = "J.s.m-6") = 31.1251e+9;
          parameter Real R_T_dorsal_pancreatic_T1_C124(unit = "J.s.m-6") = 16.6313e+9;
          parameter Real R_T_splenic_T2_C126(unit = "J.s.m-6") = 0.429028e+9;
          parameter Real R_T_hepatic_artery_proper_left_branch_C132(unit = "J.s.m-6") = 0.852085e+9;
          parameter Real R_T_hepatic_artery_proper_right_branch_C134(unit = "J.s.m-6") = 0.471438e+9;
          parameter Real R_T_middle_colic_T8_C140(unit = "J.s.m-6") = 2.6865e+9;
          parameter Real R_T_jejunal_3_T10_C144(unit = "J.s.m-6") = 1.96952e+9;
          parameter Real R_T_jejunal_6_T11_C148(unit = "J.s.m-6") = 1.96952e+9;
          parameter Real R_T_ileocolic_T9_C152(unit = "J.s.m-6") = 0.971722e+9;
          parameter Real R_T_ileal_4_T12_C156(unit = "J.s.m-6") = 1.32964e+9;
          parameter Real R_T_ileal_6_T13_C160(unit = "J.s.m-6") = 1.32964e+9;
          parameter Real R_T_superior_mesenteric_T4_C162(unit = "J.s.m-6") = 0.880745e+9;
          parameter Real R_T_inferior_segmental_T5_L170(unit = "J.s.m-6") = 0.526057e+9;
          parameter Real R_T_superior_segmental_T4_L172(unit = "J.s.m-6") = 0.526057e+9;
          parameter Real R_T_renal_posterior_branch_T3_L174(unit = "J.s.m-6") = 0.931304e+9;
          parameter Real R_T_superior_segmental_T4_R182(unit = "J.s.m-6") = 0.527823e+9;
          parameter Real R_T_inferior_segmental_T5_R184(unit = "J.s.m-6") = 0.527823e+9;
          parameter Real R_T_renal_posterior_branch_T3_R186(unit = "J.s.m-6") = 0.934431e+9;
          parameter Real R_T_inferior_mesenteric_T5_C190(unit = "J.s.m-6") = 2.16746e+9;
          parameter Real R_T_internal_iliac_T1_R218(unit = "J.s.m-6") = 0.375754e+9;
          parameter Real R_T_profundus_T2_R224(unit = "J.s.m-6") = 0.310612e+9;
          parameter Real R_T_anterior_tibial_T3_R230(unit = "J.s.m-6") = 2.22296e+9;
          parameter Real R_T_posterior_tibial_T4_R236(unit = "J.s.m-6") = 1.91761e+9;
          parameter Real R_T_internal_iliac_T1_L196(unit = "J.s.m-6") = 0.376879e+9;
          parameter Real R_T_profundus_T2_L202(unit = "J.s.m-6") = 0.310778e+9;
          parameter Real R_T_anterior_tibial_T3_L208(unit = "J.s.m-6") = 2.22365e+9;
          parameter Real R_T_posterior_tibial_T4_L214(unit = "J.s.m-6") = 1.91845e+9;
          parameter Real R_T_posterior_interosseous_T3_R40(unit = "J.s.m-6") = 4.33378e+9;
          parameter Real R_T_ulnar_T2_R42(unit = "J.s.m-6") = 1.0649e+9;
          parameter Real R_T_radial_T1_R44(unit = "J.s.m-6") = 1.04588e+9;
          parameter Real R_T_posterior_interosseous_T3_L88(unit = "J.s.m-6") = 4.34948e+9;
          parameter Real R_T_ulnar_T2_L90(unit = "J.s.m-6") = 1.08553e+9;
          parameter Real R_T_radial_T1_L92(unit = "J.s.m-6") = 1.02715e+9;
          parameter Real R_T_external_carotid_T2_R26(unit = "J.s.m-6") = 0.851183e+9;
          parameter Real R_T_external_carotid_T2_L62(unit = "J.s.m-6") = 0.854183e+9;
          parameter Real R_T_internal_carotid_R8_C(unit = "J.s.m-6") = 0.851183e+9;
          parameter Real R_T_internal_carotid_L50_C(unit = "J.s.m-6") = 0.854183e+9;
          parameter Real R_T_vertebral_R272(unit = "J.s.m-6") = 0.851183e+9;
          parameter Real R_T_vertebral_L2(unit = "J.s.m-6") = 0.854183e+9;
          parameter Real R_T_renal_R178(unit = "J.s.m-6") = 0.854183e+9;
          parameter Real R_T_renal_L166(unit = "J.s.m-6") = 0.854183e+9;
          parameter Real R_T_celiac_trunk_C116(unit = "J.s.m-6") = 0.854183e+9;
          parameter Real R_T_lateral_occipital_anterior_temporal_branch_T46_L16(unit = "J.s.m-6") = 6.17214e+9;
          parameter Real R_T_lateral_occipital_intermediate_temporal_branch_T47_L20(unit = "J.s.m-6") = 60.7317e+9;
          parameter Real R_T_lateral_occipital_posterior_temporal_branch_T48_L24(unit = "J.s.m-6") = 25.2811e+9;
          parameter Real R_T_pericallosal_parieto_occipital_branch_T60_L30(unit = "J.s.m-6") = 17.8019e+9;
          parameter Real R_T_pericallosal_precuneal_branch_T61_L34(unit = "J.s.m-6") = 12.2746e+9;
          parameter Real R_T_distal_medial_striate_T44_L40(unit = "J.s.m-6") = 65.8855e+9;
          parameter Real R_T_anterior_choroidal_T34_L116(unit = "J.s.m-6") = 9.50237e+9;
          parameter Real R_T_lateral_frontobasal_T45_L122(unit = "J.s.m-6") = 23.3019e+9;
          parameter Real R_T_prefrontal_T65_L126(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_artery_of_central_sulcus_T36_L132(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_artery_of_precentral_sulcus_T38_L134(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_polar_temporal_T63_L138(unit = "J.s.m-6") = 44.8075e+9;
          parameter Real R_T_middle_cerebral_anterior_temporal_branch_T55_L142(unit = "J.s.m-6") = 15.3827e+9;
          parameter Real R_T_middle_cerebral_middle_temporal_branch_T57_L146(unit = "J.s.m-6") = 13.5433e+9;
          parameter Real R_T_artery_of_postcentral_sulcus_T37_L150(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_middle_cerebral_posterior_temporo_occipital_branch_T59_L156(unit = "J.s.m-6") = 10.6482e+9;
          parameter Real R_T_middle_cerebral_posterior_temporal_branch_T58_L158(unit = "J.s.m-6") = 10.6482e+9;
          parameter Real R_T_anterior_parietal_T35_L162(unit = "J.s.m-6") = 15.3827e+9;
          parameter Real R_T_posterior_parietal_T64_L168(unit = "J.s.m-6") = 9.50237e+9;
          parameter Real R_T_middle_cerebral_branch_to_angular_gyrus_T56_L170(unit = "J.s.m-6") = 7.65279e+9;
          parameter Real R_T_callosomarginal_intermediomedial_frontal_branch_T41_L174(unit = "J.s.m-6") = 13.2110e+9;
          parameter Real R_T_callosomarginal_posteromedial_frontal_branch_T43_L178(unit = "J.s.m-6") = 9.39673e+9;
          parameter Real R_T_callosomarginal_cingular_branch_T40_L182(unit = "J.s.m-6") = 12.4248e+9;
          parameter Real R_T_callosomarginal_paracentral_branch_T42_L186(unit = "J.s.m-6") = 12.4248e+9;
          parameter Real R_T_medial_occipital_occipitotemporal_branch_T52_L190(unit = "J.s.m-6") = 31.0946e+9;
          parameter Real R_T_medial_occipital_parietal_branch_T53_L196(unit = "J.s.m-6") = 51.9306e+9;
          parameter Real R_T_medial_occipital_parieto_occipital_branch_T54_L198(unit = "J.s.m-6") = 10.6655e+9;
          parameter Real R_T_medial_occipital_calcarine_branch_T51_L202(unit = "J.s.m-6") = 31.0946e+9;
          parameter Real R_T_anterior_choroidal_T98_R54(unit = "J.s.m-6") = 9.50237e+9;
          parameter Real R_T_lateral_frontobasal_T110_R60(unit = "J.s.m-6") = 23.3019e+9;
          parameter Real R_T_prefrontal_T130_R64(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_artery_of_central_sulcus_T101_R70(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_artery_of_precentral_sulcus_T103_R72(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_polar_temporal_T127_R76(unit = "J.s.m-6") = 44.8075e+9;
          parameter Real R_T_middle_cerebral_anterior_temporal_branch_T119_R80(unit = "J.s.m-6") = 15.3827e+9;
          parameter Real R_T_middle_cerebral_middle_temporal_branch_T121_R84(unit = "J.s.m-6") = 13.5433e+9;
          parameter Real R_T_artery_of_postcentral_sulcus_T102_R88(unit = "J.s.m-6") = 11.9822e+9;
          parameter Real R_T_middle_cerebral_posterior_temporo_occipital_branch_T123_R94(unit = "J.s.m-6") = 10.6482e+9;
          parameter Real R_T_middle_cerebral_posterior_temporal_branch_T122_R96(unit = "J.s.m-6") = 10.6482e+9;
          parameter Real R_T_anterior_parietal_T100_R100(unit = "J.s.m-6") = 15.3827e+9;
          parameter Real R_T_posterior_parietal_T129_R106(unit = "J.s.m-6") = 9.50237e+9;
          parameter Real R_T_middle_cerebral_branch_to_angular_gyrus_T120_R108(unit = "J.s.m-6") = 7.65279e+9;
          parameter Real R_T_lateral_occipital_anterior_temporal_branch_T111_R212(unit = "J.s.m-6") = 6.17214e+9;
          parameter Real R_T_lateral_occipital_intermediate_temporal_branch_T76_R216(unit = "J.s.m-6") = 60.7317e+9;
          parameter Real R_T_lateral_occipital_posterior_temporal_branch_T112_R220(unit = "J.s.m-6") = 25.2811e+9;
          parameter Real R_T_pericallosal_parieto_occipital_branch_T124_R226(unit = "J.s.m-6") = 17.8019e+9;
          parameter Real R_T_pericallosal_precuneal_branch_T125_R230(unit = "J.s.m-6") = 12.2746e+9;
          parameter Real R_T_distal_medial_striate_T109_R236(unit = "J.s.m-6") = 65.8855e+9;
          parameter Real R_T_callosomarginal_intermediomedial_frontal_branch_T106_R242(unit = "J.s.m-6") = 13.2110e+9;
          parameter Real R_T_callosomarginal_posteromedial_frontal_branch_T108_R246(unit = "J.s.m-6") = 9.39673e+9;
          parameter Real R_T_callosomarginal_cingular_branch_T105_R250(unit = "J.s.m-6") = 12.4248e+9;
          parameter Real R_T_callosomarginal_paracentral_branch_T107_R254(unit = "J.s.m-6") = 12.4248e+9;
          parameter Real R_T_medial_occipital_occipitotemporal_branch_T80_R258(unit = "J.s.m-6") = 31.0946e+9;
          parameter Real R_T_medial_occipital_parietal_branch_T81_R264(unit = "J.s.m-6") = 51.9306e+9;
          parameter Real R_T_medial_occipital_parieto_occipital_branch_T118_R266(unit = "J.s.m-6") = 10.6655e+9;
          parameter Real R_T_medial_occipital_calcarine_branch_T79_R270(unit = "J.s.m-6") = 31.0946e+9;
          parameter Real C_T_posterior_intercostal_T1_R98(unit = "m6.J-1") = 2.2736e-12;
          parameter Real C_T_posterior_intercostal_T1_L102(unit = "m6.J-1") = 2.2162e-12;
          parameter Real C_T_posterior_intercostal_T2_R106(unit = "m6.J-1") = 2.4369e-12;
          parameter Real C_T_posterior_intercostal_T2_L110(unit = "m6.J-1") = 2.4162e-12;
          parameter Real C_T_left_gastric_T3_C120(unit = "m6.J-1") = 1.6495e-12;
          parameter Real C_T_dorsal_pancreatic_T1_C124(unit = "m6.J-1") = 15.055e-12;
          parameter Real C_T_splenic_T2_C126(unit = "m6.J-1") = 119.67e-12;
          parameter Real C_T_hepatic_artery_proper_left_branch_C132(unit = "m6.J-1") = 60.253e-12;
          parameter Real C_T_hepatic_artery_proper_right_branch_C134(unit = "m6.J-1") = 108.9e-12;
          parameter Real C_T_middle_colic_T8_C140(unit = "m6.J-1") = 19.111e-12;
          parameter Real C_T_jejunal_3_T10_C144(unit = "m6.J-1") = 26.067e-12;
          parameter Real C_T_jejunal_6_T11_C148(unit = "m6.J-1") = 26.067e-12;
          parameter Real C_T_ileocolic_T9_C152(unit = "m6.J-1") = 52.835e-12;
          parameter Real C_T_ileal_4_T12_C156(unit = "m6.J-1") = 38.612e-12;
          parameter Real C_T_ileal_6_T13_C160(unit = "m6.J-1") = 38.612e-12;
          parameter Real C_T_superior_mesenteric_T4_C162(unit = "m6.J-1") = 58.292e-12;
          parameter Real C_T_inferior_segmental_T5_L170(unit = "m6.J-1") = 97.595e-12;
          parameter Real C_T_superior_segmental_T4_L172(unit = "m6.J-1") = 97.595e-12;
          parameter Real C_T_renal_posterior_branch_T3_L174(unit = "m6.J-1") = 55.127e-12;
          parameter Real C_T_superior_segmental_T4_R182(unit = "m6.J-1") = 97.268e-12;
          parameter Real C_T_inferior_segmental_T5_R184(unit = "m6.J-1") = 97.268e-12;
          parameter Real C_T_renal_posterior_branch_T3_R186(unit = "m6.J-1") = 54.943e-12;
          parameter Real C_T_inferior_mesenteric_T5_C190(unit = "m6.J-1") = 23.687e-12;
          parameter Real C_T_internal_iliac_T1_R218(unit = "m6.J-1") = 136.63e-12;
          parameter Real C_T_profundus_T2_R224(unit = "m6.J-1") = 165.29e-12;
          parameter Real C_T_anterior_tibial_T3_R230(unit = "m6.J-1") = 23.096e-12;
          parameter Real C_T_posterior_tibial_T4_R236(unit = "m6.J-1") = 26.773e-12;
          parameter Real C_T_internal_iliac_T1_L196(unit = "m6.J-1") = 136.23e-12;
          parameter Real C_T_profundus_T2_L202(unit = "m6.J-1") = 165.2e-12;
          parameter Real C_T_anterior_tibial_T3_L208(unit = "m6.J-1") = 23.088e-12;
          parameter Real C_T_posterior_tibial_T4_L214(unit = "m6.J-1") = 26.761e-12;
          parameter Real C_T_posterior_interosseous_T3_R40(unit = "m6.J-1") = 11.847e-12;
          parameter Real C_T_ulnar_T2_R42(unit = "m6.J-1") = 48.212e-12;
          parameter Real C_T_radial_T1_R44(unit = "m6.J-1") = 49.088e-12;
          parameter Real C_T_posterior_interosseous_T3_L88(unit = "m6.J-1") = 11.804e-12;
          parameter Real C_T_ulnar_T2_L90(unit = "m6.J-1") = 47.295e-12;
          parameter Real C_T_radial_T1_L92(unit = "m6.J-1") = 49.983e-12;
          parameter Real C_T_external_carotid_T2_R26(unit = "m6.J-1") = 60.317e-12;
          parameter Real C_T_external_carotid_T2_L62(unit = "m6.J-1") = 60.105e-12;
          parameter Real C_T_internal_carotid_R8_C(unit = "m6.J-1") = 60.317e-12;
          parameter Real C_T_internal_carotid_L50_C(unit = "m6.J-1") = 60.105e-12;
          parameter Real C_T_vertebral_R272(unit = "m6.J-1") = 60.317e-12;
          parameter Real C_T_vertebral_L2(unit = "m6.J-1") = 60.105e-12;
          parameter Real C_T_renal_R178(unit = "m6.J-1") = 60.105e-12;
          parameter Real C_T_renal_L166(unit = "m6.J-1") = 60.105e-12;
          parameter Real C_T_celiac_trunk_C116(unit = "m6.J-1") = 60.105e-12;
          parameter Real C_T_lateral_occipital_anterior_temporal_branch_T46_L16(unit = "m6.J-1") = 2.60095e-12;
          parameter Real C_T_lateral_occipital_intermediate_temporal_branch_T47_L20(unit = "m6.J-1") = 0.264334e-12;
          parameter Real C_T_lateral_occipital_posterior_temporal_branch_T48_L24(unit = "m6.J-1") = 0.634997e-12;
          parameter Real C_T_pericallosal_parieto_occipital_branch_T60_L30(unit = "m6.J-1") = 0.901782e-12;
          parameter Real C_T_pericallosal_precuneal_branch_T61_L34(unit = "m6.J-1") = 1.30786e-12;
          parameter Real C_T_distal_medial_striate_T44_L40(unit = "m6.J-1") = 0.243656e-12;
          parameter Real C_T_anterior_choroidal_T34_L116(unit = "m6.J-1") = 1.68941e-12;
          parameter Real C_T_lateral_frontobasal_T45_L122(unit = "m6.J-1") = 0.688932e-12;
          parameter Real C_T_prefrontal_T65_L126(unit = "m6.J-1") = 1.33978e-12;
          parameter Real C_T_artery_of_central_sulcus_T36_L132(unit = "m6.J-1") = 1.33978e-12;
          parameter Real C_T_artery_of_precentral_sulcus_T38_L134(unit = "m6.J-1") = 1.33978e-12;
          parameter Real C_T_polar_temporal_T63_L138(unit = "m6.J-1") = 0.35828e-12;
          parameter Real C_T_middle_cerebral_anterior_temporal_branch_T55_L142(unit = "m6.J-1") = 1.0436e-12;
          parameter Real C_T_middle_cerebral_middle_temporal_branch_T57_L146(unit = "m6.J-1") = 1.1853e-12;
          parameter Real C_T_artery_of_postcentral_sulcus_T37_L150(unit = "m6.J-1") = 1.3398e-12;
          parameter Real C_T_middle_cerebral_posterior_temporo_occipital_branch_T59_L156(unit = "m6.J-1") = 1.5076e-12;
          parameter Real C_T_middle_cerebral_posterior_temporal_branch_T58_L158(unit = "m6.J-1") = 1.5076e-12;
          parameter Real C_T_anterior_parietal_T35_L162(unit = "m6.J-1") = 1.0436e-12;
          parameter Real C_T_posterior_parietal_T64_L168(unit = "m6.J-1") = 1.6894e-12;
          parameter Real C_T_middle_cerebral_branch_to_angular_gyrus_T56_L170(unit = "m6.J-1") = 2.0977e-12;
          parameter Real C_T_callosomarginal_intermediomedial_frontal_branch_T41_L174(unit = "m6.J-1") = 1.2152e-12;
          parameter Real C_T_callosomarginal_posteromedial_frontal_branch_T43_L178(unit = "m6.J-1") = 1.7084e-12;
          parameter Real C_T_callosomarginal_cingular_branch_T40_L182(unit = "m6.J-1") = 1.2921e-12;
          parameter Real C_T_callosomarginal_paracentral_branch_T42_L186(unit = "m6.J-1") = 1.2921e-12;
          parameter Real C_T_medial_occipital_occipitotemporal_branch_T52_L190(unit = "m6.J-1") = 0.51628e-12;
          parameter Real C_T_medial_occipital_parietal_branch_T53_L196(unit = "m6.J-1") = 0.30913e-12;
          parameter Real C_T_medial_occipital_parieto_occipital_branch_T54_L198(unit = "m6.J-1") = 1.5052e-12;
          parameter Real C_T_medial_occipital_calcarine_branch_T51_L202(unit = "m6.J-1") = 0.51628e-12;
          parameter Real C_T_anterior_choroidal_T98_R54(unit = "m6.J-1") = 1.6894e-12;
          parameter Real C_T_lateral_frontobasal_T110_R60(unit = "m6.J-1") = 0.68893e-12;
          parameter Real C_T_prefrontal_T130_R64(unit = "m6.J-1") = 1.3398e-12;
          parameter Real C_T_artery_of_central_sulcus_T101_R70(unit = "m6.J-1") = 1.3398e-12;
          parameter Real C_T_artery_of_precentral_sulcus_T103_R72(unit = "m6.J-1") = 1.3398e-12;
          parameter Real C_T_polar_temporal_T127_R76(unit = "m6.J-1") = 0.35828e-12;
          parameter Real C_T_middle_cerebral_anterior_temporal_branch_T119_R80(unit = "m6.J-1") = 1.0436e-12;
          parameter Real C_T_middle_cerebral_middle_temporal_branch_T121_R84(unit = "m6.J-1") = 1.1853e-12;
          parameter Real C_T_artery_of_postcentral_sulcus_T102_R88(unit = "m6.J-1") = 1.3398e-12;
          parameter Real C_T_middle_cerebral_posterior_temporo_occipital_branch_T123_R94(unit = "m6.J-1") = 1.5076e-12;
          parameter Real C_T_middle_cerebral_posterior_temporal_branch_T122_R96(unit = "m6.J-1") = 1.5076e-12;
          parameter Real C_T_anterior_parietal_T100_R100(unit = "m6.J-1") = 1.0436e-12;
          parameter Real C_T_posterior_parietal_T129_R106(unit = "m6.J-1") = 1.6894e-12;
          parameter Real C_T_middle_cerebral_branch_to_angular_gyrus_T120_R108(unit = "m6.J-1") = 2.0977e-12;
          parameter Real C_T_lateral_occipital_anterior_temporal_branch_T111_R212(unit = "m6.J-1") = 2.6009e-12;
          parameter Real C_T_lateral_occipital_intermediate_temporal_branch_T76_R216(unit = "m6.J-1") = 0.26433e-12;
          parameter Real C_T_lateral_occipital_posterior_temporal_branch_T112_R220(unit = "m6.J-1") = 0.635e-12;
          parameter Real C_T_pericallosal_parieto_occipital_branch_T124_R226(unit = "m6.J-1") = 0.90178e-12;
          parameter Real C_T_pericallosal_precuneal_branch_T125_R230(unit = "m6.J-1") = 1.3079e-12;
          parameter Real C_T_distal_medial_striate_T109_R236(unit = "m6.J-1") = 0.24366e-12;
          parameter Real C_T_callosomarginal_intermediomedial_frontal_branch_T106_R242(unit = "m6.J-1") = 1.2152e-12;
          parameter Real C_T_callosomarginal_posteromedial_frontal_branch_T108_R246(unit = "m6.J-1") = 1.7084e-12;
          parameter Real C_T_callosomarginal_cingular_branch_T105_R250(unit = "m6.J-1") = 1.2921e-12;
          parameter Real C_T_callosomarginal_paracentral_branch_T107_R254(unit = "m6.J-1") = 1.2921e-12;
          parameter Real C_T_medial_occipital_occipitotemporal_branch_T80_R258(unit = "m6.J-1") = 0.51628e-12;
          parameter Real C_T_medial_occipital_parietal_branch_T81_R264(unit = "m6.J-1") = 0.30913e-12;
          parameter Real C_T_medial_occipital_parieto_occipital_branch_T118_R266(unit = "m6.J-1") = 1.5052e-12;
          parameter Real C_T_medial_occipital_calcarine_branch_T79_R270(unit = "m6.J-1") = 0.51628e-12;
          parameter Real E_ascending_aorta_A(unit = "Pa") = 0.4e+6;
          parameter Real E_ascending_aorta_B(unit = "Pa") = 0.4e+6;
          parameter Real E_ascending_aorta_C(unit = "Pa") = 0.4e+6;
          parameter Real E_ascending_aorta_D(unit = "Pa") = 0.4e+6;
          parameter Real E_aortic_arch_C2(unit = "Pa") = 0.4e+6;
          parameter Real E_brachiocephalic_trunk_C4(unit = "Pa") = 0.4e+6;
          parameter Real E_aortic_arch_C46(unit = "Pa") = 0.4e+6;
          parameter Real E_aortic_arch_C64(unit = "Pa") = 0.4e+6;
          parameter Real E_aortic_arch_C94(unit = "Pa") = 0.4e+6;
          parameter Real E_thoracic_aorta_C96(unit = "Pa") = 0.4e+6;
          parameter Real E_thoracic_aorta_C100(unit = "Pa") = 0.4e+6;
          parameter Real E_thoracic_aorta_C104(unit = "Pa") = 0.4e+6;
          parameter Real E_thoracic_aorta_C108(unit = "Pa") = 0.4e+6;
          parameter Real E_thoracic_aorta_C112(unit = "Pa") = 0.4e+6;
          parameter Real E_abdominal_aorta_C114(unit = "Pa") = 0.4e+6;
          parameter Real E_abdominal_aorta_C136(unit = "Pa") = 0.4e+6;
          parameter Real E_abdominal_aorta_C164(unit = "Pa") = 0.4e+6;
          parameter Real E_abdominal_aorta_C176(unit = "Pa") = 0.4e+6;
          parameter Real E_abdominal_aorta_C188(unit = "Pa") = 0.4e+6;
          parameter Real E_abdominal_aorta_C192(unit = "Pa") = 0.4e+6;
          parameter Real E_posterior_intercostal_T1_R98(unit = "Pa") = 0.4e+6;
          parameter Real E_posterior_intercostal_T1_L102(unit = "Pa") = 0.4e+6;
          parameter Real E_posterior_intercostal_T2_R106(unit = "Pa") = 0.4e+6;
          parameter Real E_posterior_intercostal_T2_L110(unit = "Pa") = 0.4e+6;
          parameter Real E_celiac_trunk_C116(unit = "Pa") = 0.4e+6;
          parameter Real E_splenic_T2_C118(unit = "Pa") = 0.4e+6;
          parameter Real E_left_gastric_T3_C120(unit = "Pa") = 0.4e+6;
          parameter Real E_splenic_T2_C122(unit = "Pa") = 0.4e+6;
          parameter Real E_dorsal_pancreatic_T1_C124(unit = "Pa") = 0.4e+6;
          parameter Real E_splenic_T2_C126(unit = "Pa") = 0.4e+6;
          parameter Real E_common_hepatic_C128(unit = "Pa") = 0.4e+6;
          parameter Real E_hepatic_artery_proper_C130(unit = "Pa") = 0.4e+6;
          parameter Real E_hepatic_artery_proper_left_branch_C132(unit = "Pa") = 0.4e+6;
          parameter Real E_hepatic_artery_proper_right_branch_C134(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_mesenteric_T4_C138(unit = "Pa") = 0.4e+6;
          parameter Real E_middle_colic_T8_C140(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_mesenteric_T4_C142(unit = "Pa") = 0.4e+6;
          parameter Real E_jejunal_3_T10_C144(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_mesenteric_T4_C146(unit = "Pa") = 0.4e+6;
          parameter Real E_jejunal_6_T11_C148(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_mesenteric_T4_C150(unit = "Pa") = 0.4e+6;
          parameter Real E_ileocolic_T9_C152(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_mesenteric_T4_C154(unit = "Pa") = 0.4e+6;
          parameter Real E_ileal_4_T12_C156(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_mesenteric_T4_C158(unit = "Pa") = 0.4e+6;
          parameter Real E_ileal_6_T13_C160(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_mesenteric_T4_C162(unit = "Pa") = 0.4e+6;
          parameter Real E_renal_L166(unit = "Pa") = 0.4e+6;
          parameter Real E_renal_anterior_branch_L168(unit = "Pa") = 0.4e+6;
          parameter Real E_inferior_segmental_T5_L170(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_segmental_T4_L172(unit = "Pa") = 0.4e+6;
          parameter Real E_renal_posterior_branch_T3_L174(unit = "Pa") = 0.4e+6;
          parameter Real E_renal_R178(unit = "Pa") = 0.4e+6;
          parameter Real E_renal_anterior_branch_R180(unit = "Pa") = 0.4e+6;
          parameter Real E_superior_segmental_T4_R182(unit = "Pa") = 0.4e+6;
          parameter Real E_inferior_segmental_T5_R184(unit = "Pa") = 0.4e+6;
          parameter Real E_renal_posterior_branch_T3_R186(unit = "Pa") = 0.4e+6;
          parameter Real E_inferior_mesenteric_T5_C190(unit = "Pa") = 0.4e+6;
          parameter Real E_common_iliac_R216(unit = "Pa") = 0.4e+6;
          parameter Real E_internal_iliac_T1_R218(unit = "Pa") = 1.6e+6;
          parameter Real E_external_iliac_R220(unit = "Pa") = 0.8e+6;
          parameter Real E_femoral_R222(unit = "Pa") = 0.8e+6;
          parameter Real E_profundus_T2_R224(unit = "Pa") = 0.8e+6;
          parameter Real E_femoral_R226(unit = "Pa") = 0.8e+6;
          parameter Real E_popliteal_R228(unit = "Pa") = 0.8e+6;
          parameter Real E_anterior_tibial_T3_R230(unit = "Pa") = 1.6e+6;
          parameter Real E_popliteal_R232(unit = "Pa") = 1.6e+6;
          parameter Real E_tibiofibular_trunk_R234(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_tibial_T4_R236(unit = "Pa") = 1.6e+6;
          parameter Real E_common_iliac_L194(unit = "Pa") = 0.4e+6;
          parameter Real E_internal_iliac_T1_L196(unit = "Pa") = 1.6e+6;
          parameter Real E_external_iliac_L198(unit = "Pa") = 0.8e+6;
          parameter Real E_femoral_L200(unit = "Pa") = 0.8e+6;
          parameter Real E_profundus_T2_L202(unit = "Pa") = 0.8e+6;
          parameter Real E_femoral_L204(unit = "Pa") = 0.8e+6;
          parameter Real E_popliteal_L206(unit = "Pa") = 0.8e+6;
          parameter Real E_anterior_tibial_T3_L208(unit = "Pa") = 1.6e+6;
          parameter Real E_popliteal_L210(unit = "Pa") = 1.6e+6;
          parameter Real E_tibiofibular_trunk_L212(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_tibial_T4_L214(unit = "Pa") = 1.6e+6;
          parameter Real E_subclavian_R28(unit = "Pa") = 0.4e+6;
          parameter Real E_subclavian_R30(unit = "Pa") = 0.4e+6;
          parameter Real E_axillary_R32(unit = "Pa") = 0.4e+6;
          parameter Real E_brachial_R34(unit = "Pa") = 0.4e+6;
          parameter Real E_ulnar_T2_R36(unit = "Pa") = 0.8e+6;
          parameter Real E_common_interosseous_R38(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_interosseous_T3_R40(unit = "Pa") = 1.6e+6;
          parameter Real E_ulnar_T2_R42(unit = "Pa") = 0.8e+6;
          parameter Real E_radial_T1_R44(unit = "Pa") = 0.8e+6;
          parameter Real E_subclavian_L66(unit = "Pa") = 0.4e+6;
          parameter Real E_subclavian_L78(unit = "Pa") = 0.4e+6;
          parameter Real E_axillary_L80(unit = "Pa") = 0.4e+6;
          parameter Real E_brachial_L82(unit = "Pa") = 0.4e+6;
          parameter Real E_ulnar_T2_L84(unit = "Pa") = 0.8e+6;
          parameter Real E_common_interosseous_L86(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_interosseous_T3_L88(unit = "Pa") = 1.6e+6;
          parameter Real E_ulnar_T2_L90(unit = "Pa") = 0.8e+6;
          parameter Real E_radial_T1_L92(unit = "Pa") = 0.8e+6;
          parameter Real E_common_carotid_R6_A(unit = "Pa") = 0.2e+6;
          parameter Real E_common_carotid_R6_B(unit = "Pa") = 0.2e+6;
          parameter Real E_common_carotid_R6_C(unit = "Pa") = 0.2e+6;
          parameter Real E_internal_carotid_R8_A(unit = "Pa") = 0.8e+6;
          parameter Real E_internal_carotid_R8_B(unit = "Pa") = 0.8e+6;
          parameter Real E_internal_carotid_R8_C(unit = "Pa") = 1.6e+6;
          parameter Real E_external_carotid_T2_R26(unit = "Pa") = 0.8e+6;
          parameter Real E_common_carotid_L48_A(unit = "Pa") = 0.2e+6;
          parameter Real E_common_carotid_L48_B(unit = "Pa") = 0.2e+6;
          parameter Real E_common_carotid_L48_C(unit = "Pa") = 0.2e+6;
          parameter Real E_common_carotid_L48_D(unit = "Pa") = 0.2e+6;
          parameter Real E_internal_carotid_L50_A(unit = "Pa") = 0.8e+6;
          parameter Real E_internal_carotid_L50_B(unit = "Pa") = 0.8e+6;
          parameter Real E_internal_carotid_L50_C(unit = "Pa") = 1.6e+6;
          parameter Real E_external_carotid_T2_L62(unit = "Pa") = 0.8e+6;
          parameter Real E_vertebral_L2(unit = "Pa") = 0.8e+6;
          parameter Real E_vertebral_R272(unit = "Pa") = 0.8e+6;
          parameter Real E_basilar_C4(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_cerebral_precommunicating_part_L6(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_cerebral_precommunicating_part_R204(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_communicating_L8(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_communicating_R206(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_cerebral_postcommunicating_part_L12(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_cerebral_postcommunicating_part_R208(unit = "Pa") = 1.6e+6;
          parameter Real E_occipital_lateral_L14(unit = "Pa") = 1.6e+6;
          parameter Real E_occipital_lateral_R210(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_occipital_anterior_temporal_branch_T111_R212(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_occipital_anterior_temporal_branch_T46_L16(unit = "Pa") = 1.6e+6;
          parameter Real E_occipital_lateral_L18(unit = "Pa") = 1.6e+6;
          parameter Real E_occipital_lateral_R214(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_occipital_intermediate_temporal_branch_T47_L20(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_occipital_intermediate_temporal_branch_T76_R216(unit = "Pa") = 1.6e+6;
          parameter Real E_occipital_lateral_L22(unit = "Pa") = 1.6e+6;
          parameter Real E_occipital_lateral_R218(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_occipital_posterior_temporal_branch_T48_L24(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_occipital_posterior_temporal_branch_T112_R220(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_L26(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_R222(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_dorsal_branch_to_corpus_callosum_L28(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_dorsal_branch_to_corpus_callosum_R224(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_parieto_occipital_branch_T60_L30(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_parieto_occipital_branch_T124_R226(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_L32(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_R228(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_precuneal_branch_T61_L34(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_precuneal_branch_T125_R230(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_L36(unit = "Pa") = 1.6e+6;
          parameter Real E_pericallosal_R232(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_cerebral_L38(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_cerebral_R234(unit = "Pa") = 1.6e+6;
          parameter Real E_distal_medial_striate_T44_L40(unit = "Pa") = 1.6e+6;
          parameter Real E_distal_medial_striate_T109_R236(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_cerebral_L42(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_cerebral_R238(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_communicating_C44(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_cerebral_L110(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_cerebral_R46(unit = "Pa") = 1.6e+6;
          parameter Real E_internal_carotid_L112(unit = "Pa") = 1.6e+6;
          parameter Real E_internal_carotid_R48(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_L114(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_R52(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_choroidal_T34_L116(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_choroidal_T98_R54(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_L118(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_R56(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_superior_terminal_branch_L120(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_superior_terminal_branch_R58(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_frontobasal_T45_L122(unit = "Pa") = 1.6e+6;
          parameter Real E_lateral_frontobasal_T110_R60(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_superior_terminal_branch_L124(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_superior_terminal_branch_R62(unit = "Pa") = 1.6e+6;
          parameter Real E_prefrontal_T65_L126(unit = "Pa") = 1.6e+6;
          parameter Real E_prefrontal_T130_R64(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_superior_terminal_branch_L128(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_superior_terminal_branch_R66(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_precentral_sulcus_T38_L130(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_precentral_sulcus_T103_R68(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_central_sulcus_T36_L132(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_central_sulcus_T101_R70(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_precentral_sulcus_T38_L134(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_precentral_sulcus_T103_R72(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_L136(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_R74(unit = "Pa") = 1.6e+6;
          parameter Real E_polar_temporal_T63_L138(unit = "Pa") = 1.6e+6;
          parameter Real E_polar_temporal_T127_R76(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_L140(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_R78(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_anterior_temporal_branch_T55_L142(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_anterior_temporal_branch_T119_R80(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_L144(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_R82(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_middle_temporal_branch_T57_L146(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_middle_temporal_branch_T121_R84(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_L148(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_R86(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_postcentral_sulcus_T37_L150(unit = "Pa") = 1.6e+6;
          parameter Real E_artery_of_postcentral_sulcus_T102_R88(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_L152(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_R90(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_posterior_temporal_branch_T58_L154(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_posterior_temporal_branch_T122_R92(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_posterior_temporo_occipital_branch_T59_L156(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_posterior_temporo_occipital_branch_T123_R94(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_posterior_temporal_branch_T58_L158(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_posterior_temporal_branch_T122_R96(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_L160(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_R98(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_parietal_T35_L162(unit = "Pa") = 1.6e+6;
          parameter Real E_anterior_parietal_T100_R100(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_L164(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_inferior_terminal_branch_R102(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_branch_to_angular_gyrus_T56_L166(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_branch_to_angular_gyrus_T120_R104(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_parietal_T64_L168(unit = "Pa") = 1.6e+6;
          parameter Real E_posterior_parietal_T129_R106(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_branch_to_angular_gyrus_T56_L170(unit = "Pa") = 1.6e+6;
          parameter Real E_middle_cerebral_branch_to_angular_gyrus_T120_R108(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_L172(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_R240(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_intermediomedial_frontal_branch_T41_L174(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_intermediomedial_frontal_branch_T106_R242(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_L176(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_R244(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_posteromedial_frontal_branch_T43_L178(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_posteromedial_frontal_branch_T108_R246(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_L180(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_R248(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_cingular_branch_T40_L182(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_cingular_branch_T105_R250(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_L184(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_R252(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_paracentral_branch_T42_L186(unit = "Pa") = 1.6e+6;
          parameter Real E_callosomarginal_paracentral_branch_T107_R254(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_L188(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_R256(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_occipitotemporal_branch_T52_L190(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_occipitotemporal_branch_T80_R258(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_L192(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_R260(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_parieto_occipital_branch_T54_L194(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_parieto_occipital_branch_T118_R262(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_parietal_branch_T53_L196(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_parietal_branch_T81_R264(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_parieto_occipital_branch_T54_L198(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_parieto_occipital_branch_T118_R266(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_L200(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_R268(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_calcarine_branch_T51_L202(unit = "Pa") = 1.6e+6;
          parameter Real E_medial_occipital_calcarine_branch_T79_R270(unit = "Pa") = 1.6e+6;
        equation

        end Parameters_Systemic;

        model Parameters_Venous
          parameter Real r_superior_vena_cava_C2(unit = "m") = 0.975e-2;
          parameter Real r_azygos_vein_T1_C4(unit = "m") = 0.38e-2;
          parameter Real r_superior_vena_cava_C88(unit = "m") = 0.975e-2;
          parameter Real r_azygos_vein_T1_C86(unit = "m") = 0.38e-2;
          parameter Real r_inferior_vena_cava_C8(unit = "m") = 0.995e-2;
          parameter Real r_hepatic_vein_T1_C10(unit = "m") = 0.6965e-2;
          parameter Real r_inferior_vena_cava_C12(unit = "m") = 0.98728e-2;
          parameter Real r_venous_perforator_T2_C14(unit = "m") = 0.4975e-2;
          parameter Real r_inferior_vena_cava_C16(unit = "m") = 0.94787e-2;
          parameter Real r_renal_vein_T1_R18(unit = "m") = 0.6e-2;
          parameter Real r_inferior_vena_cava_C20(unit = "m") = 0.94787e-2;
          parameter Real r_renal_vein_T1_L22(unit = "m") = 0.55e-2;
          parameter Real r_inferior_vena_cava_C24(unit = "m") = 0.94661e-2;
          parameter Real r_common_iliac_vein_L56(unit = "m") = 0.6314e-2;
          parameter Real r_common_iliac_vein_R26(unit = "m") = 0.6314e-2;
          parameter Real r_external_iliac_vein_R28(unit = "m") = 0.63e-2;
          parameter Real r_internal_iliac_vein_T1_R30(unit = "m") = 0.491e-2;
          parameter Real r_external_iliac_vein_R32(unit = "m") = 0.63e-2;
          parameter Real r_femoral_vein_R34(unit = "m") = 0.568e-2;
          parameter Real r_great_saphenous_vein_T7_R36(unit = "m") = 0.181e-2;
          parameter Real r_femoral_vein_R38(unit = "m") = 0.568e-2;
          parameter Real r_profunda_femoris_vein_T2_R40(unit = "m") = 0.4e-2;
          parameter Real r_femoral_vein_R42(unit = "m") = 0.568e-2;
          parameter Real r_venous_perforator_T3_R44(unit = "m") = 0.2e-2;
          parameter Real r_femoral_vein_R46(unit = "m") = 0.568e-2;
          parameter Real r_popliteal_vein_R48(unit = "m") = 0.381e-2;
          parameter Real r_anterior_tibial_vein_T4_R50(unit = "m") = 0.1175e-2;
          parameter Real r_popliteal_vein_R52(unit = "m") = 0.381e-2;
          parameter Real r_posterior_tibial_vein_T6_R54(unit = "m") = 0.17e-2;
          parameter Real r_external_iliac_vein_L58(unit = "m") = 0.63e-2;
          parameter Real r_internal_iliac_vein_T1_L60(unit = "m") = 0.491e-2;
          parameter Real r_external_iliac_vein_L62(unit = "m") = 0.63e-2;
          parameter Real r_femoral_vein_L64(unit = "m") = 0.568e-2;
          parameter Real r_great_saphenous_vein_T7_L66(unit = "m") = 0.181e-2;
          parameter Real r_femoral_vein_L68(unit = "m") = 0.568e-2;
          parameter Real r_profunda_femoris_vein_T2_L70(unit = "m") = 0.4e-2;
          parameter Real r_femoral_vein_L72(unit = "m") = 0.568e-2;
          parameter Real r_venous_perforator_T3_L74(unit = "m") = 0.568e-2;
          parameter Real r_femoral_vein_L76(unit = "m") = 0.568e-2;
          parameter Real r_popliteal_vein_L78(unit = "m") = 0.381e-2;
          parameter Real r_anterior_tibial_vein_T4_L80(unit = "m") = 0.1175e-2;
          parameter Real r_popliteal_vein_L82(unit = "m") = 0.381e-2;
          parameter Real r_posterior_tibial_vein_T6_L84(unit = "m") = 0.17e-2;
          parameter Real r_brachiocephalic_vein_R90(unit = "m") = 0.8e-2;
          parameter Real r_brachiocephalic_vein_L124(unit = "m") = 0.75e-2;
          parameter Real r_vertebral_vein_R92(unit = "m") = 0.27e-2;
          parameter Real r_brachiocephalic_vein_R94(unit = "m") = 0.8e-2;
          parameter Real r_subclavian_vein_R96(unit = "m") = 0.56e-2;
          parameter Real r_internal_jugular_vein_R122(unit = "m") = 0.75e-2;
          parameter Real r_external_jugular_vein_R98(unit = "m") = 0.225e-2;
          parameter Real r_subclavian_vein_R100(unit = "m") = 0.56e-2;
          parameter Real r_axillary_vein_R102(unit = "m") = 0.56e-2;
          parameter Real r_brachial_vein_R104(unit = "m") = 0.185e-2;
          parameter Real r_brachial_vein_R114(unit = "m") = 0.17e-2;
          parameter Real r_venous_perforator_T2_R106(unit = "m") = 0.17e-2;
          parameter Real r_brachial_vein_R108(unit = "m") = 0.185e-2;
          parameter Real r_ulnar_vein_T7_R110(unit = "m") = 0.085e-2;
          parameter Real r_brachial_vein_R112(unit = "m") = 0.185e-2;
          parameter Real r_venous_perforator_T1_R116(unit = "m") = 0.185e-2;
          parameter Real r_brachial_vein_R118(unit = "m") = 0.17e-2;
          parameter Real r_radial_vein_T3_R120(unit = "m") = 0.0665e-2;
          parameter Real r_vertebral_vein_L126(unit = "m") = 0.27e-2;
          parameter Real r_brachiocephalic_vein_L128(unit = "m") = 0.75e-2;
          parameter Real r_subclavian_vein_L130(unit = "m") = 0.56e-2;
          parameter Real r_internal_jugular_vein_L156(unit = "m") = 0.6e-2;
          parameter Real r_external_jugular_vein_L132(unit = "m") = 0.225e-2;
          parameter Real r_subclavian_vein_L134(unit = "m") = 0.56e-2;
          parameter Real r_axillary_vein_L136(unit = "m") = 0.56e-2;
          parameter Real r_brachial_vein_L138(unit = "m") = 0.185e-2;
          parameter Real r_brachial_vein_L148(unit = "m") = 0.17e-2;
          parameter Real r_venous_perforator_T2_L140(unit = "m") = 0.185e-2;
          parameter Real r_brachial_vein_L142(unit = "m") = 0.185e-2;
          parameter Real r_ulnar_vein_T7_L144(unit = "m") = 0.085e-2;
          parameter Real r_brachial_vein_L146(unit = "m") = 0.185e-2;
          parameter Real r_venous_perforator_T1_L150(unit = "m") = 0.17e-2;
          parameter Real r_brachial_vein_L152(unit = "m") = 0.17e-2;
          parameter Real r_radial_vein_T3_L154(unit = "m") = 0.0665e-2;
          parameter Real h_superior_vena_cava_C2(unit = "m") = 0.11778e-2;
          parameter Real h_azygos_vein_T1_C4(unit = "m") = 0.06384e-2;
          parameter Real h_superior_vena_cava_C88(unit = "m") = 0.11778e-2;
          parameter Real h_azygos_vein_T1_C86(unit = "m") = 0.06384e-2;
          parameter Real h_inferior_vena_cava_C8(unit = "m") = 0.11974e-2;
          parameter Real h_hepatic_vein_T1_C10(unit = "m") = 0.09111e-2;
          parameter Real h_inferior_vena_cava_C12(unit = "m") = 0.11899e-2;
          parameter Real h_venous_perforator_T2_C14(unit = "m") = 0.0736e-2;
          parameter Real h_inferior_vena_cava_C16(unit = "m") = 0.11513e-2;
          parameter Real h_renal_vein_T1_R18(unit = "m") = 0.08241e-2;
          parameter Real h_inferior_vena_cava_C20(unit = "m") = 0.11513e-2;
          parameter Real h_renal_vein_T1_L22(unit = "m") = 0.07806e-2;
          parameter Real h_inferior_vena_cava_C24(unit = "m") = 0.11501e-2;
          parameter Real h_common_iliac_vein_L56(unit = "m") = 0.0852e-2;
          parameter Real h_common_iliac_vein_R26(unit = "m") = 0.0852e-2;
          parameter Real h_external_iliac_vein_R28(unit = "m") = 0.08508e-2;
          parameter Real h_internal_iliac_vein_T1_R30(unit = "m") = 0.07306e-2;
          parameter Real h_external_iliac_vein_R32(unit = "m") = 0.08508e-2;
          parameter Real h_femoral_vein_R34(unit = "m") = 0.07962e-2;
          parameter Real h_great_saphenous_vein_T7_R36(unit = "m") = 0.04381e-2;
          parameter Real h_femoral_vein_R38(unit = "m") = 0.07962e-2;
          parameter Real h_profunda_femoris_vein_T2_R40(unit = "m") = 0.0655e-2;
          parameter Real h_femoral_vein_R42(unit = "m") = 0.07962e-2;
          parameter Real h_venous_perforator_T3_R44(unit = "m") = 0.0463e-2;
          parameter Real h_femoral_vein_R46(unit = "m") = 0.07962e-2;
          parameter Real h_popliteal_vein_R48(unit = "m") = 0.06392e-2;
          parameter Real h_anterior_tibial_vein_T4_R50(unit = "m") = 0.03354e-2;
          parameter Real h_popliteal_vein_R52(unit = "m") = 0.06392e-2;
          parameter Real h_posterior_tibial_vein_T6_R54(unit = "m") = 0.04226e-2;
          parameter Real h_external_iliac_vein_L58(unit = "m") = 0.08508e-2;
          parameter Real h_internal_iliac_vein_T1_L60(unit = "m") = 0.07306e-2;
          parameter Real h_external_iliac_vein_L62(unit = "m") = 0.08508e-2;
          parameter Real h_femoral_vein_L64(unit = "m") = 0.07962e-2;
          parameter Real h_great_saphenous_vein_T7_L66(unit = "m") = 0.04381e-2;
          parameter Real h_femoral_vein_L68(unit = "m") = 0.07962e-2;
          parameter Real h_profunda_femoris_vein_T2_L70(unit = "m") = 0.0655e-2;
          parameter Real h_femoral_vein_L72(unit = "m") = 0.07962e-2;
          parameter Real h_venous_perforator_T3_L74(unit = "m") = 0.07962e-2;
          parameter Real h_femoral_vein_L76(unit = "m") = 0.07962e-2;
          parameter Real h_popliteal_vein_L78(unit = "m") = 0.06392e-2;
          parameter Real h_anterior_tibial_vein_T4_L80(unit = "m") = 0.03354e-2;
          parameter Real h_popliteal_vein_L82(unit = "m") = 0.06392e-2;
          parameter Real h_posterior_tibial_vein_T6_L84(unit = "m") = 0.04226e-2;
          parameter Real h_brachiocephalic_vein_R90(unit = "m") = 0.10082e-2;
          parameter Real h_brachiocephalic_vein_L124(unit = "m") = 0.09609e-2;
          parameter Real h_vertebral_vein_R92(unit = "m") = 0.05402e-2;
          parameter Real h_brachiocephalic_vein_R94(unit = "m") = 0.10082e-2;
          parameter Real h_subclavian_vein_R96(unit = "m") = 0.07892e-2;
          parameter Real h_internal_jugular_vein_R122(unit = "m") = 0.09609e-2;
          parameter Real h_external_jugular_vein_R98(unit = "m") = 0.04928e-2;
          parameter Real h_subclavian_vein_R100(unit = "m") = 0.07892e-2;
          parameter Real h_axillary_vein_R102(unit = "m") = 0.07892e-2;
          parameter Real h_brachial_vein_R104(unit = "m") = 0.04435e-2;
          parameter Real h_brachial_vein_R114(unit = "m") = 0.04226e-2;
          parameter Real h_venous_perforator_T2_R106(unit = "m") = 0.04226e-2;
          parameter Real h_brachial_vein_R108(unit = "m") = 0.04435e-2;
          parameter Real h_ulnar_vein_T7_R110(unit = "m") = 0.02665e-2;
          parameter Real h_brachial_vein_R112(unit = "m") = 0.04435e-2;
          parameter Real h_venous_perforator_T1_R116(unit = "m") = 0.04435e-2;
          parameter Real h_brachial_vein_R118(unit = "m") = 0.04226e-2;
          parameter Real h_radial_vein_T3_R120(unit = "m") = 0.02206e-2;
          parameter Real h_vertebral_vein_L126(unit = "m") = 0.05402e-2;
          parameter Real h_brachiocephalic_vein_L128(unit = "m") = 0.09609e-2;
          parameter Real h_subclavian_vein_L130(unit = "m") = 0.07892e-2;
          parameter Real h_internal_jugular_vein_L156(unit = "m") = 0.08241e-2;
          parameter Real h_external_jugular_vein_L132(unit = "m") = 0.04928e-2;
          parameter Real h_subclavian_vein_L134(unit = "m") = 0.07892e-2;
          parameter Real h_axillary_vein_L136(unit = "m") = 0.07892e-2;
          parameter Real h_brachial_vein_L138(unit = "m") = 0.04435e-2;
          parameter Real h_brachial_vein_L148(unit = "m") = 0.04226e-2;
          parameter Real h_venous_perforator_T2_L140(unit = "m") = 0.04435e-2;
          parameter Real h_brachial_vein_L142(unit = "m") = 0.04435e-2;
          parameter Real h_ulnar_vein_T7_L144(unit = "m") = 0.02665e-2;
          parameter Real h_brachial_vein_L146(unit = "m") = 0.04435e-2;
          parameter Real h_venous_perforator_T1_L150(unit = "m") = 0.04226e-2;
          parameter Real h_brachial_vein_L152(unit = "m") = 0.04226e-2;
          parameter Real h_radial_vein_T3_L154(unit = "m") = 0.02206e-2;
          parameter Real l_superior_vena_cava_C2(unit = "m") = 1.19864e-2;
          parameter Real l_azygos_vein_T1_C4(unit = "m") = 24.1235e-2;
          parameter Real l_superior_vena_cava_C88(unit = "m") = 2.06549e-2;
          parameter Real l_azygos_vein_T1_C86(unit = "m") = 4.16471e-2;
          parameter Real l_inferior_vena_cava_C8(unit = "m") = 1.72606e-2;
          parameter Real l_hepatic_vein_T1_C10(unit = "m") = 1.38104e-2;
          parameter Real l_inferior_vena_cava_C12(unit = "m") = 4.75594e-2;
          parameter Real l_venous_perforator_T2_C14(unit = "m") = 1.61605e-2;
          parameter Real l_inferior_vena_cava_C16(unit = "m") = 4.05694e-2;
          parameter Real l_renal_vein_T1_R18(unit = "m") = 2.98282e-2;
          parameter Real l_inferior_vena_cava_C20(unit = "m") = 0.283038e-2;
          parameter Real l_renal_vein_T1_L22(unit = "m") = 3.13626e-2;
          parameter Real l_inferior_vena_cava_C24(unit = "m") = 11.3168e-2;
          parameter Real l_common_iliac_vein_L56(unit = "m") = 6.24014e-2;
          parameter Real l_common_iliac_vein_R26(unit = "m") = 5.73581e-2;
          parameter Real l_external_iliac_vein_R28(unit = "m") = 0.840236e-2;
          parameter Real l_internal_iliac_vein_T1_R30(unit = "m") = 5.91456e-2;
          parameter Real l_external_iliac_vein_R32(unit = "m") = 9.74495e-2;
          parameter Real l_femoral_vein_R34(unit = "m") = 0.521142e-2;
          parameter Real l_great_saphenous_vein_T7_R36(unit = "m") = 92.6589e-2;
          parameter Real l_femoral_vein_R38(unit = "m") = 4.13055e-2;
          parameter Real l_profunda_femoris_vein_T2_R40(unit = "m") = 22.9655e-2;
          parameter Real l_femoral_vein_R42(unit = "m") = 29.2094e-2;
          parameter Real l_venous_perforator_T3_R44(unit = "m") = 1.01988e-2;
          parameter Real l_femoral_vein_R46(unit = "m") = 1.60864e-2;
          parameter Real l_popliteal_vein_R48(unit = "m") = 14.5387e-2;
          parameter Real l_anterior_tibial_vein_T4_R50(unit = "m") = 38.09e-2;
          parameter Real l_popliteal_vein_R52(unit = "m") = 3.2338e-2;
          parameter Real l_posterior_tibial_vein_T6_R54(unit = "m") = 30.5697e-2;
          parameter Real l_external_iliac_vein_L58(unit = "m") = 0.841264e-2;
          parameter Real l_internal_iliac_vein_T1_L60(unit = "m") = 6.08336e-2;
          parameter Real l_external_iliac_vein_L62(unit = "m") = 9.56715e-2;
          parameter Real l_femoral_vein_L64(unit = "m") = 0.609248e-2;
          parameter Real l_great_saphenous_vein_T7_L66(unit = "m") = 92.6588e-2;
          parameter Real l_femoral_vein_L68(unit = "m") = 4.13055e-2;
          parameter Real l_profunda_femoris_vein_T2_L70(unit = "m") = 22.9654e-2;
          parameter Real l_femoral_vein_L72(unit = "m") = 29.2094e-2;
          parameter Real l_venous_perforator_T3_L74(unit = "m") = 1.01987e-2;
          parameter Real l_femoral_vein_L76(unit = "m") = 1.60863e-2;
          parameter Real l_popliteal_vein_L78(unit = "m") = 14.5387e-2;
          parameter Real l_anterior_tibial_vein_T4_L80(unit = "m") = 38.09e-2;
          parameter Real l_popliteal_vein_L82(unit = "m") = 3.23382e-2;
          parameter Real l_posterior_tibial_vein_T6_L84(unit = "m") = 30.5698e-2;
          parameter Real l_brachiocephalic_vein_R90(unit = "m") = 3.60711e-2;
          parameter Real l_brachiocephalic_vein_L124(unit = "m") = 7.69083e-2;
          parameter Real l_vertebral_vein_R92(unit = "m") = 19.7909e-2;
          parameter Real l_brachiocephalic_vein_R94(unit = "m") = 0.817729e-2;
          parameter Real l_subclavian_vein_R96(unit = "m") = 0.811432e-2;
          parameter Real l_internal_jugular_vein_R122(unit = "m") = 17.8797e-2;
          parameter Real l_external_jugular_vein_R98(unit = "m") = 13.807e-2;
          parameter Real l_subclavian_vein_R100(unit = "m") = 3.31851e-2;
          parameter Real l_axillary_vein_R102(unit = "m") = 11.6597e-2;
          parameter Real l_brachial_vein_R104(unit = "m") = 20.4344e-2;
          parameter Real l_brachial_vein_R114(unit = "m") = 21.6418e-2;
          parameter Real l_venous_perforator_T2_R106(unit = "m") = 1.64071e-2;
          parameter Real l_brachial_vein_R108(unit = "m") = 2.81322e-2;
          parameter Real l_ulnar_vein_T7_R110(unit = "m") = 26.5376e-2;
          parameter Real l_brachial_vein_R112(unit = "m") = 1.06711e-2;
          parameter Real l_venous_perforator_T1_R116(unit = "m") = 1.56083e-2;
          parameter Real l_brachial_vein_R118(unit = "m") = 2.085e-2;
          parameter Real l_radial_vein_T3_R120(unit = "m") = 25.8545e-2;
          parameter Real l_vertebral_vein_L126(unit = "m") = 18.5123e-2;
          parameter Real l_brachiocephalic_vein_L128(unit = "m") = 0.470158e-2;
          parameter Real l_subclavian_vein_L130(unit = "m") = 0.721557e-2;
          parameter Real l_internal_jugular_vein_L156(unit = "m") = 16.9554e-2;
          parameter Real l_external_jugular_vein_L132(unit = "m") = 13.5402e-2;
          parameter Real l_subclavian_vein_L134(unit = "m") = 3.23403e-2;
          parameter Real l_axillary_vein_L136(unit = "m") = 11.7916e-2;
          parameter Real l_brachial_vein_L138(unit = "m") = 20.4344e-2;
          parameter Real l_brachial_vein_L148(unit = "m") = 21.6418e-2;
          parameter Real l_venous_perforator_T2_L140(unit = "m") = 1.64073e-2;
          parameter Real l_brachial_vein_L142(unit = "m") = 2.81325e-2;
          parameter Real l_ulnar_vein_T7_L144(unit = "m") = 26.5376e-2;
          parameter Real l_brachial_vein_L146(unit = "m") = 1.06711e-2;
          parameter Real l_venous_perforator_T1_L150(unit = "m") = 1.56083e-2;
          parameter Real l_brachial_vein_L152(unit = "m") = 2.085e-2;
          parameter Real l_radial_vein_T3_L154(unit = "m") = 25.8545e-2;
          parameter Real C_superior_vena_cava_C2(unit = "m6.J-1") = 0.555603e-10;
          parameter Real C_azygos_vein_T1_C4(unit = "m6.J-1") = 1.22148e-10;
          parameter Real C_superior_vena_cava_C88(unit = "m6.J-1") = 0.957416e-10;
          parameter Real C_azygos_vein_T1_C86(unit = "m6.J-1") = 0.210878e-10;
          parameter Real C_inferior_vena_cava_C8(unit = "m6.J-1") = 0.829339e-10;
          parameter Real C_hepatic_vein_T1_C10(unit = "m6.J-1") = 0.30168e-10;
          parameter Real C_inferior_vena_cava_C12(unit = "m6.J-1") = 2.21275e-10;
          parameter Real C_venous_perforator_T2_C14(unit = "m6.J-1") = 0.159254e-10;
          parameter Real C_inferior_vena_cava_C16(unit = "m6.J-1") = 1.80509e-10;
          parameter Real C_renal_vein_T1_R18(unit = "m6.J-1") = 0.460512e-10;
          parameter Real C_inferior_vena_cava_C20(unit = "m6.J-1") = 0.123145e-10;
          parameter Real C_renal_vein_T1_L22(unit = "m6.J-1") = 0.393746e-10;
          parameter Real C_inferior_vena_cava_C24(unit = "m6.J-1") = 4.63354e-10;
          parameter Real C_common_iliac_vein_L56(unit = "m6.J-1") = 1.08597e-10;
          parameter Real C_common_iliac_vein_R26(unit = "m6.J-1") = 0.998205e-10;
          parameter Real C_external_iliac_vein_R28(unit = "m6.J-1") = 0.14547e-10;
          parameter Real C_internal_iliac_vein_T1_R30(unit = "m6.J-1") = 0.564484e-10;
          parameter Real C_external_iliac_vein_R32(unit = "m6.J-1") = 1.68715e-10;
          parameter Real C_femoral_vein_R34(unit = "m6.J-1") = 0.0706573e-10;
          parameter Real C_great_saphenous_vein_T7_R36(unit = "m6.J-1") = 0.738807e-10;
          parameter Real C_femoral_vein_R38(unit = "m6.J-1") = 0.56e-10;
          parameter Real C_profunda_femoris_vein_T2_R40(unit = "m6.J-1") = 1.32176e-10;
          parameter Real C_femoral_vein_R42(unit = "m6.J-1") = 3.96026e-10;
          parameter Real C_venous_perforator_T3_R44(unit = "m6.J-1") = 0.00550536e-10;
          parameter Real C_femoral_vein_R46(unit = "m6.J-1") = 0.218102e-10;
          parameter Real C_popliteal_vein_R48(unit = "m6.J-1") = 0.741013e-10;
          parameter Real C_anterior_tibial_vein_T4_R50(unit = "m6.J-1") = 0.108531e-10;
          parameter Real C_popliteal_vein_R52(unit = "m6.J-1") = 0.164821e-10;
          parameter Real C_posterior_tibial_vein_T6_R54(unit = "m6.J-1") = 0.209333e-10;
          parameter Real C_external_iliac_vein_L58(unit = "m6.J-1") = 0.145648e-10;
          parameter Real C_internal_iliac_vein_T1_L60(unit = "m6.J-1") = 0.580595e-10;
          parameter Real C_external_iliac_vein_L62(unit = "m6.J-1") = 1.65636e-10;
          parameter Real C_femoral_vein_L64(unit = "m6.J-1") = 0.0826029e-10;
          parameter Real C_great_saphenous_vein_T7_L66(unit = "m6.J-1") = 0.738807e-10;
          parameter Real C_femoral_vein_L68(unit = "m6.J-1") = 0.560026e-10;
          parameter Real C_profunda_femoris_vein_T2_L70(unit = "m6.J-1") = 1.32176e-10;
          parameter Real C_femoral_vein_L72(unit = "m6.J-1") = 3.96026e-10;
          parameter Real C_venous_perforator_T3_L74(unit = "m6.J-1") = 0.138275e-10;
          parameter Real C_femoral_vein_L76(unit = "m6.J-1") = 0.218101e-10;
          parameter Real C_popliteal_vein_L78(unit = "m6.J-1") = 0.741013e-10;
          parameter Real C_anterior_tibial_vein_T4_L80(unit = "m6.J-1") = 0.108531e-10;
          parameter Real C_popliteal_vein_L82(unit = "m6.J-1") = 0.164822e-10;
          parameter Real C_posterior_tibial_vein_T6_L84(unit = "m6.J-1") = 0.209333e-10;
          parameter Real C_brachiocephalic_vein_R90(unit = "m6.J-1") = 1.07899e-10;
          parameter Real C_brachiocephalic_vein_L124(unit = "m6.J-1") = 1.98898e-10;
          parameter Real C_vertebral_vein_R92(unit = "m6.J-1") = 0.424745e-10;
          parameter Real C_brachiocephalic_vein_R94(unit = "m6.J-1") = 0.244606e-10;
          parameter Real C_subclavian_vein_R96(unit = "m6.J-1") = 0.106357e-10;
          parameter Real C_internal_jugular_vein_R122(unit = "m6.J-1") = 3.36401e-10;
          parameter Real C_external_jugular_vein_R98(unit = "m6.J-1") = 0.187997e-10;
          parameter Real C_subclavian_vein_R100(unit = "m6.J-1") = 0.434967e-10;
          parameter Real C_axillary_vein_R102(unit = "m6.J-1") = 1.52827e-10;
          parameter Real C_brachial_vein_R104(unit = "m6.J-1") = 0.171851e-10;
          parameter Real C_brachial_vein_R114(unit = "m6.J-1") = 0.148197e-10;
          parameter Real C_venous_perforator_T2_R106(unit = "m6.J-1") = 0.0112351e-10;
          parameter Real C_brachial_vein_R108(unit = "m6.J-1") = 0.0236588e-10;
          parameter Real C_ulnar_vein_T7_R110(unit = "m6.J-1") = 0.0360248e-10;
          parameter Real C_brachial_vein_R112(unit = "m6.J-1") = 0.00897427e-10;
          parameter Real C_venous_perforator_T1_R116(unit = "m6.J-1") = 0.0131264e-10;
          parameter Real C_brachial_vein_R118(unit = "m6.J-1") = 0.0142775e-10;
          parameter Real C_radial_vein_T3_R120(unit = "m6.J-1") = 0.0203056e-10;
          parameter Real C_vertebral_vein_L126(unit = "m6.J-1") = 0.397304e-10;
          parameter Real C_brachiocephalic_vein_L128(unit = "m6.J-1") = 0.121591e-10;
          parameter Real C_subclavian_vein_L130(unit = "m6.J-1") = 0.0945766e-10;
          parameter Real C_internal_jugular_vein_L156(unit = "m6.J-1") = 1.92288e-10;
          parameter Real C_external_jugular_vein_L132(unit = "m6.J-1") = 0.184364e-10;
          parameter Real C_subclavian_vein_L134(unit = "m6.J-1") = 0.423894e-10;
          parameter Real C_axillary_vein_L136(unit = "m6.J-1") = 1.54556e-10;
          parameter Real C_brachial_vein_L138(unit = "m6.J-1") = 0.171851e-10;
          parameter Real C_brachial_vein_L148(unit = "m6.J-1") = 0.148197e-10;
          parameter Real C_venous_perforator_T2_L140(unit = "m6.J-1") = 0.0137984e-10;
          parameter Real C_brachial_vein_L142(unit = "m6.J-1") = 0.0236591e-10;
          parameter Real C_ulnar_vein_T7_L144(unit = "m6.J-1") = 0.0360249e-10;
          parameter Real C_brachial_vein_L146(unit = "m6.J-1") = 0.00897427e-10;
          parameter Real C_venous_perforator_T1_L150(unit = "m6.J-1") = 0.0106881e-10;
          parameter Real C_brachial_vein_L152(unit = "m6.J-1") = 0.0142775e-10;
          parameter Real C_radial_vein_T3_L154(unit = "m6.J-1") = 0.0203056e-10;
          parameter Real E_superior_vena_cava_C2(unit = "Pa") = 0.555603e+6;
          parameter Real E_azygos_vein_T1_C4(unit = "Pa") = 1.22148e+6;
          parameter Real E_superior_vena_cava_C88(unit = "Pa") = 0.957416e+6;
          parameter Real E_azygos_vein_T1_C86(unit = "Pa") = 0.210878e+6;
          parameter Real E_inferior_vena_cava_C8(unit = "Pa") = 0.829339e+6;
          parameter Real E_hepatic_vein_T1_C10(unit = "Pa") = 0.30168e+6;
          parameter Real E_inferior_vena_cava_C12(unit = "Pa") = 2.21275e+6;
          parameter Real E_venous_perforator_T2_C14(unit = "Pa") = 0.159254e+6;
          parameter Real E_inferior_vena_cava_C16(unit = "Pa") = 1.80509e+6;
          parameter Real E_renal_vein_T1_R18(unit = "Pa") = 0.460512e+6;
          parameter Real E_inferior_vena_cava_C20(unit = "Pa") = 0.123145e+6;
          parameter Real E_renal_vein_T1_L22(unit = "Pa") = 0.393746e+6;
          parameter Real E_inferior_vena_cava_C24(unit = "Pa") = 4.63354e+6;
          parameter Real E_common_iliac_vein_L56(unit = "Pa") = 1.08597e+6;
          parameter Real E_common_iliac_vein_R26(unit = "Pa") = 0.998205e+6;
          parameter Real E_external_iliac_vein_R28(unit = "Pa") = 0.14547e+6;
          parameter Real E_internal_iliac_vein_T1_R30(unit = "Pa") = 0.564484e+6;
          parameter Real E_external_iliac_vein_R32(unit = "Pa") = 1.68715e+6;
          parameter Real E_femoral_vein_R34(unit = "Pa") = 0.0706573e+6;
          parameter Real E_great_saphenous_vein_T7_R36(unit = "Pa") = 0.738807e+6;
          parameter Real E_femoral_vein_R38(unit = "Pa") = 0.56e+6;
          parameter Real E_profunda_femoris_vein_T2_R40(unit = "Pa") = 1.32176e+6;
          parameter Real E_femoral_vein_R42(unit = "Pa") = 3.96026e+6;
          parameter Real E_venous_perforator_T3_R44(unit = "Pa") = 0.00550536e+6;
          parameter Real E_femoral_vein_R46(unit = "Pa") = 0.218102e+6;
          parameter Real E_popliteal_vein_R48(unit = "Pa") = 0.741013e+6;
          parameter Real E_anterior_tibial_vein_T4_R50(unit = "Pa") = 0.108531e+6;
          parameter Real E_popliteal_vein_R52(unit = "Pa") = 0.164821e+6;
          parameter Real E_posterior_tibial_vein_T6_R54(unit = "Pa") = 0.209333e+6;
          parameter Real E_external_iliac_vein_L58(unit = "Pa") = 0.145648e+6;
          parameter Real E_internal_iliac_vein_T1_L60(unit = "Pa") = 0.580595e+6;
          parameter Real E_external_iliac_vein_L62(unit = "Pa") = 1.65636e+6;
          parameter Real E_femoral_vein_L64(unit = "Pa") = 0.0826029e+6;
          parameter Real E_great_saphenous_vein_T7_L66(unit = "Pa") = 0.738807e+6;
          parameter Real E_femoral_vein_L68(unit = "Pa") = 0.560026e+6;
          parameter Real E_profunda_femoris_vein_T2_L70(unit = "Pa") = 1.32176e+6;
          parameter Real E_femoral_vein_L72(unit = "Pa") = 3.96026e+6;
          parameter Real E_venous_perforator_T3_L74(unit = "Pa") = 0.138275e+6;
          parameter Real E_femoral_vein_L76(unit = "Pa") = 0.218101e+6;
          parameter Real E_popliteal_vein_L78(unit = "Pa") = 0.741013e+6;
          parameter Real E_anterior_tibial_vein_T4_L80(unit = "Pa") = 0.108531e+6;
          parameter Real E_popliteal_vein_L82(unit = "Pa") = 0.164822e+6;
          parameter Real E_posterior_tibial_vein_T6_L84(unit = "Pa") = 0.209333e+6;
          parameter Real E_brachiocephalic_vein_R90(unit = "Pa") = 1.07899e+6;
          parameter Real E_brachiocephalic_vein_L124(unit = "Pa") = 1.98898e+6;
          parameter Real E_vertebral_vein_R92(unit = "Pa") = 0.424745e+6;
          parameter Real E_brachiocephalic_vein_R94(unit = "Pa") = 0.244606e+6;
          parameter Real E_subclavian_vein_R96(unit = "Pa") = 0.106357e+6;
          parameter Real E_internal_jugular_vein_R122(unit = "Pa") = 3.36401e+6;
          parameter Real E_external_jugular_vein_R98(unit = "Pa") = 0.187997e+6;
          parameter Real E_subclavian_vein_R100(unit = "Pa") = 0.434967e+6;
          parameter Real E_axillary_vein_R102(unit = "Pa") = 1.52827e+6;
          parameter Real E_brachial_vein_R104(unit = "Pa") = 0.171851e+6;
          parameter Real E_brachial_vein_R114(unit = "Pa") = 0.148197e+6;
          parameter Real E_venous_perforator_T2_R106(unit = "Pa") = 0.0112351e+6;
          parameter Real E_brachial_vein_R108(unit = "Pa") = 0.0236588e+6;
          parameter Real E_ulnar_vein_T7_R110(unit = "Pa") = 0.0360248e+6;
          parameter Real E_brachial_vein_R112(unit = "Pa") = 0.00897427e+6;
          parameter Real E_venous_perforator_T1_R116(unit = "Pa") = 0.0131264e+6;
          parameter Real E_brachial_vein_R118(unit = "Pa") = 0.0142775e+6;
          parameter Real E_radial_vein_T3_R120(unit = "Pa") = 0.0203056e+6;
          parameter Real E_vertebral_vein_L126(unit = "Pa") = 0.397304e+6;
          parameter Real E_brachiocephalic_vein_L128(unit = "Pa") = 0.121591e+6;
          parameter Real E_subclavian_vein_L130(unit = "Pa") = 0.0945766e+6;
          parameter Real E_internal_jugular_vein_L156(unit = "Pa") = 1.92288e+6;
          parameter Real E_external_jugular_vein_L132(unit = "Pa") = 0.184364e+6;
          parameter Real E_subclavian_vein_L134(unit = "Pa") = 0.423894e+6;
          parameter Real E_axillary_vein_L136(unit = "Pa") = 1.54556e+6;
          parameter Real E_brachial_vein_L138(unit = "Pa") = 0.171851e+6;
          parameter Real E_brachial_vein_L148(unit = "Pa") = 0.148197e+6;
          parameter Real E_venous_perforator_T2_L140(unit = "Pa") = 0.0137984e+6;
          parameter Real E_brachial_vein_L142(unit = "Pa") = 0.0236591e+6;
          parameter Real E_ulnar_vein_T7_L144(unit = "Pa") = 0.0360249e+6;
          parameter Real E_brachial_vein_L146(unit = "Pa") = 0.00897427e+6;
          parameter Real E_venous_perforator_T1_L150(unit = "Pa") = 0.0106881e+6;
          parameter Real E_brachial_vein_L152(unit = "Pa") = 0.0142775e+6;
          parameter Real E_radial_vein_T3_L154(unit = "Pa") = 0.0203056e+6;
        equation

        end Parameters_Venous;
      end Parameters_Venous_cellml;

        model Pulmonary
          ADAN_main.Components.AdanVenousRed._b580e.Parameters_Venous_cellml.Parameters_Pulmonary
            Parameters_Pulmonary1 annotation (Placement(transformation(extent={
                    {-100,80},{-80,100}})));
          parameter Boolean UseSimplifiedInertance = true;
          input Real t(unit = "s");
          Real C_pas(unit = "m6.J-1");
          Real C_pat(unit = "m6.J-1");
          Real C_par(unit = "m6.J-1");
          Real C_pcp(unit = "m6.J-1");
          Real C_pvn(unit = "m6.J-1");
          Real C_pvc(unit = "m6.J-1");
          Real R_pas(unit = "J.s.m-6");
          Real R_pat(unit = "J.s.m-6");
          Real R_par(unit = "J.s.m-6");
          Real R_pcp(unit = "J.s.m-6");
          Real R_psh(unit = "J.s.m-6");
          Real R_pvn(unit = "J.s.m-6");
          Real I_pas(unit = "J.s2.m-6");
          Real I_pat(unit = "J.s2.m-6");
          Real I_par(unit = "J.s2.m-6");
          Real I_pcp(unit = "J.s2.m-6");
          Real I_pvn(unit = "J.s2.m-6");
          Real I_psh(unit = "J.s2.m-6");
          Real u_pas(unit = "Pa", start = 4000.0);
          Real u_pat(unit = "Pa", start = 0.0);
          Real u_par(unit = "Pa", start = 0.0);
          Real u_pcp(unit = "Pa", start = 0.0);
          Real u_pvn(unit = "Pa", start = 0.0);
          input Real u_la(unit = "Pa");
          Real v_pas(unit = "m3.s-1", nominal = 1e-6);
          Real v_pat(unit = "m3.s-1", nominal = 1e-6);
          Real v_par(unit = "m3.s-1", nominal = 1e-6);
          Real v_pcp(unit = "m3.s-1", nominal = 1e-6);
          Real v_psh(unit = "m3.s-1", nominal = 1e-6);
          Real v_pvn(unit = "m3.s-1", nominal = 1e-6);
          input Real v_puv(unit = "m3.s-1");
          Physiolibrary.Types.Volume V_pas = (u_pas - thoracic_pressure)*C_pas "Pulmonary artery volume";
          Physiolibrary.Types.Volume V_pat = (u_pat - thoracic_pressure)*C_pat "Pulmonary arteries volume";
          Physiolibrary.Types.Volume V_par = (u_par - thoracic_pressure)*C_par "Pulmonary arterioles volume";
          Physiolibrary.Types.Volume V_pcp = (u_pcp - thoracic_pressure)*C_pcp "Pulmonary capillaries volume";
          Physiolibrary.Types.Volume V_pvn = (u_pvn - thoracic_pressure)*C_pvn "Pulmonary veins volume";
          Physiolibrary.Types.Volume total_stressed_volume = V_pas + V_pat + V_par + V_pcp + V_pvn;
          input Physiolibrary.Types.Pressure thoracic_pressure;
        //  Physiolibrary.Types.Volume V_pvc = u_pvc*C_pvc;
        equation
          C_pas = Parameters_Pulmonary1.C_pas;
          C_pat = Parameters_Pulmonary1.C_pat;
          C_par = Parameters_Pulmonary1.C_par;
          C_pcp = Parameters_Pulmonary1.C_pcp;
          C_pvn = Parameters_Pulmonary1.C_pvn;
          C_pvc = Parameters_Pulmonary1.C_pvc;
          R_pas = Parameters_Pulmonary1.R_pas;
          R_pat = Parameters_Pulmonary1.R_pat;
          R_par = Parameters_Pulmonary1.R_par;
          R_pcp = Parameters_Pulmonary1.R_pcp;
          R_psh = Parameters_Pulmonary1.R_psh;
          R_pvn = Parameters_Pulmonary1.R_pvn;
          I_pas = Parameters_Pulmonary1.I_pas;
          I_pat = Parameters_Pulmonary1.I_pat;
          I_par = Parameters_Pulmonary1.I_par;
          I_pcp = Parameters_Pulmonary1.I_pcp;
          I_pvn = Parameters_Pulmonary1.I_pvn;
          I_psh = Parameters_Pulmonary1.I_psh;

              der(V_pas) = (v_puv-v_pas);
              der(V_pat) = (v_pas-v_pat);
              der(V_par) = (v_pat-v_psh-v_par);
              der(V_pcp) = (v_par-v_pcp);
              der(V_pvn) = (v_pcp+v_psh-v_pvn);
              der(v_pas) = (u_pas-u_pat-v_pas*R_pas)/I_pas;
              der(v_pat) = (u_pat-u_par-v_pat*R_pat)/I_pat;

              if not UseSimplifiedInertance then
                der(v_par) = (u_par-u_pcp-v_par*R_par)/I_par;
                der(v_pcp) = (u_pcp-u_pvn-v_pcp*R_pcp)/I_pcp;
                der(v_pvn) = (u_pvn-u_la-v_pvn*R_pvn)/I_pvn;
                der(v_psh) = (u_par-u_pvn-v_psh*R_psh)/I_psh;
              else
                0 = (u_par-u_pcp-v_par*R_par);
                0 = (u_pcp-u_pvn-v_pcp*R_pcp);
                0 = (u_pvn-u_la-v_pvn*R_pvn);
                0 = (u_par-u_pvn-v_psh*R_psh);
              end if;


        end Pulmonary;
      end _b580e;

      package _7af7a4

      package Parameters86_cellml

        model Parameters_Heart
          parameter Real T(unit = "s") = 1.0;
          parameter Real CQ_trv = 34.6427e-6 "UnitValve";
          parameter Real CQ_puv = 30.3124e-6 "UnitValve";
          parameter Real CQ_miv = 34.6427e-6 "UnitValve";
          parameter Real CQ_aov = 30.3124e-6 "UnitValve";
          parameter Real E_la_max(unit = "J.m-6") = 213.28e+5;
          parameter Real E_lv_max(unit = "J.m-6") = 2857.39e+5;
          parameter Real E_ra_max(unit = "J.m-6") = 173.29e+5;
          parameter Real E_rv_max(unit = "J.m-6") = 799.8e+5;
          parameter Real E_la_min(unit = "J.m-6") = 119.97e+5;
          parameter Real E_lv_min(unit = "J.m-6") = 106.64e+5;
          parameter Real E_ra_min(unit = "J.m-6") = 93.31e+5;
          parameter Real E_rv_min(unit = "J.m-6") = 66.65e+5;
          parameter Real q_lv_0(unit = "m3") = 5.0e-6;
          parameter Real q_rv_0(unit = "m3") = 10.0e-6;
          parameter Real q_la_0(unit = "m3") = 4.0e-6;
          parameter Real q_ra_0(unit = "m3") = 4.0e-6;
          parameter Real tau_1_lv(unit = "1") = 0.215;
          parameter Real tau_1_rv(unit = "1") = 0.215;
          parameter Real tau_1_la(unit = "1") = 0.042;
          parameter Real tau_1_ra(unit = "1") = 0.042;
          parameter Real tau_2_lv(unit = "1") = 0.362;
          parameter Real tau_2_rv(unit = "1") = 0.362;
          parameter Real tau_2_la(unit = "1") = 0.138;
          parameter Real tau_2_ra(unit = "1") = 0.138;
          parameter Real t_onset_a(unit = "1") = 0.65;
          parameter Real m_1_lv(unit = "1") = 1.32;
          parameter Real m_1_rv(unit = "1") = 1.32;
          parameter Real m_1_la(unit = "1") = 1.99;
          parameter Real m_1_ra(unit = "1") = 1.99;
          parameter Real m_2_lv(unit = "1") = 21.9;
          parameter Real m_2_rv(unit = "1") = 21.9;
          parameter Real m_2_la(unit = "1") = 11.2;
          parameter Real m_2_ra(unit = "1") = 11.2;
        equation

        end Parameters_Heart;
      end Parameters86_cellml;

        model HeartInputs
          ADAN_main.Components.AdanVenousRed._7af7a4.Parameters86_cellml.Parameters_Heart
            Parameters_Heart1
            annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
          input Real t(unit = "s");
          input Real thoracic_pressure;
          Real mt(unit = "s", start = 0, fixed = true);
        //   Real mt_(unit = "s");
          Real mta(unit = "s", start = 0, fixed = true);
        //   Real mta_;
          input Real T(unit = "s");
          Physiolibrary.Types.Frequency f = 1/T;
          discrete Modelica.SIunits.Time last_beat(start = 0);
        //  Real int_f;
          Real d_mt = t_onset_a/f;

          Real CQ_trv;//(unit = "UnitValve");
          Real CQ_puv;//(unit = "UnitValve");
          Real CQ_miv;//(unit = "UnitValve");
          Real CQ_aov;//(unit = "UnitValve");
          Real q_ra_0(unit = "m3");
          Real q_rv_0(unit = "m3");
          Real q_la_0(unit = "m3");
          Real q_lv_0(unit = "m3");
          Real E_lv_max(unit = "J.m-6");
          Real E_lv_min(unit = "J.m-6");
          Real E_la_max(unit = "J.m-6");
          Real E_la_min(unit = "J.m-6");
          Real E_rv_max(unit = "J.m-6");
          Real E_rv_min(unit = "J.m-6");
          Real E_ra_max(unit = "J.m-6");
          Real E_ra_min(unit = "J.m-6");
          Real tau_1_lv(unit = "1");
          Real tau_2_lv(unit = "1");
          Real tau_1_la(unit = "1");
          Real tau_2_la(unit = "1");
          Real tau_1_rv(unit = "1");
          Real tau_2_rv(unit = "1");
          Real tau_1_ra(unit = "1");
          Real tau_2_ra(unit = "1");
          Real m_1_lv(unit = "1");
          Real m_2_lv(unit = "1");
          Real m_1_la(unit = "1");
          Real m_2_la(unit = "1");
          Real m_1_rv(unit = "1");
          Real m_2_rv(unit = "1");
          Real m_1_ra(unit = "1");
          Real m_2_ra(unit = "1");
          Real g_1_lv(unit = "1");
          Real g_2_lv(unit = "1");
          Real g_1_la(unit = "1");
          Real g_2_la(unit = "1");
          Real g_1_rv(unit = "1");
          Real g_2_rv(unit = "1");
          Real g_1_ra(unit = "1");
          Real g_2_ra(unit = "1");
          Real t_onset_a(unit = "1") = 0.85 "According to Mynard thesis 2011 and in contrdiction to Safaei repos";
          Real E_lv(unit = "J.m-6");
          Real E_la(unit = "J.m-6");
          Real E_rv(unit = "J.m-6");
          Real E_ra(unit = "J.m-6");
          Real H_la(unit = "1");
          Real H_lv(unit = "1");
          Real H_ra(unit = "1");
          Real H_rv(unit = "1");
          Real u_ra(unit = "Pa");
          Real u_rv(unit = "Pa");
          Real u_la(unit = "Pa");
          Real u_lv(unit = "Pa");
          input Real u_sas(unit = "Pa");
          input Real u_par(unit = "Pa");
          Real v_trv(unit = "m3.s-1");
          Real v_puv(unit = "m3.s-1");
          Real v_miv(unit = "m3.s-1");
          Real v_aov(unit = "m3.s-1");
          input Real v_sup_venacava(unit = "m3.s-1");
          input Real v_inf_venacava(unit = "m3.s-1");
          input Real v_pvn(unit = "m3.s-1");
          Real q_lv(unit = "m3", start = 600.0e-6);
          Real q_rv(unit = "m3", start = 600.0e-6);
          Real q_la(unit = "m3", start = 20.0e-6);
          Real q_ra(unit = "m3", start = 20.0e-6);
          Physiolibrary.Types.Volume stressed_volume = q_la + q_lv + q_ra + q_rv;
          Physiolibrary.Types.Volume zp_volume = q_la_0 + q_lv_0 + q_ra_0 + q_rv_0;
          Physiolibrary.Types.Volume total_volume = stressed_volume + zp_volume;
          Physiolibrary.Types.VolumeFlowRate v_ra_filling = v_sup_venacava+v_inf_venacava;
          Physiolibrary.Types.Pressure boundedTP_ra = max(thoracic_pressure*min(q_ra/q_ra_0, 1), 0);
          Physiolibrary.Types.Pressure boundedTP_rv = max(thoracic_pressure*min(q_rv/q_rv_0, 1), 0);
          Physiolibrary.Types.Pressure boundedTP_la = max(thoracic_pressure*min(q_la/q_la_0, 1), 0);
          Physiolibrary.Types.Pressure boundedTP_lv = max(thoracic_pressure*min(q_lv/q_lv_0, 1), 0);
        equation
        //  T = Parameters_Heart1.T;
          CQ_trv = Parameters_Heart1.CQ_trv;
          CQ_puv = Parameters_Heart1.CQ_puv;
          CQ_miv = Parameters_Heart1.CQ_miv;
          CQ_aov = Parameters_Heart1.CQ_aov;
          E_lv_max = Parameters_Heart1.E_lv_max;
          E_lv_min = Parameters_Heart1.E_lv_min;
          E_la_max = Parameters_Heart1.E_la_max;
          E_la_min = Parameters_Heart1.E_la_min;
          E_rv_max = Parameters_Heart1.E_rv_max;
          E_rv_min = Parameters_Heart1.E_rv_min;
          E_ra_max = Parameters_Heart1.E_ra_max;
          E_ra_min = Parameters_Heart1.E_ra_min;
          tau_1_lv = Parameters_Heart1.tau_1_lv;
          tau_2_lv = Parameters_Heart1.tau_2_lv;
          tau_1_la = Parameters_Heart1.tau_1_la;
          tau_2_la = Parameters_Heart1.tau_2_la;
          tau_1_rv = Parameters_Heart1.tau_1_rv;
          tau_2_rv = Parameters_Heart1.tau_2_rv;
          tau_1_ra = Parameters_Heart1.tau_1_ra;
          tau_2_ra = Parameters_Heart1.tau_2_ra;
          m_1_lv = Parameters_Heart1.m_1_lv;
          m_2_lv = Parameters_Heart1.m_2_lv;
          m_1_la = Parameters_Heart1.m_1_la;
          m_2_la = Parameters_Heart1.m_2_la;
          m_1_rv = Parameters_Heart1.m_1_rv;
          m_2_rv = Parameters_Heart1.m_2_rv;
          m_1_ra = Parameters_Heart1.m_1_ra;
          m_2_ra = Parameters_Heart1.m_2_ra;
        //  t_onset_a = Parameters_Heart1.t_onset_a;
          q_ra_0 = Parameters_Heart1.q_ra_0;
          q_rv_0 = Parameters_Heart1.q_rv_0;
          q_la_0 = Parameters_Heart1.q_la_0;
          q_lv_0 = Parameters_Heart1.q_lv_0;

              H_lv = pow((tau_1_lv+tau_2_lv)/(tau_1_lv+tau_1_lv), m_1_lv)/(pow((tau_1_lv+tau_2_lv)/(tau_1_lv+tau_1_lv), m_1_lv)+1)*1/(pow((tau_1_lv+tau_2_lv)/(tau_2_lv+tau_2_lv), m_2_lv)+1);
              H_la = pow((tau_1_la+tau_2_la)/(tau_1_la+tau_1_la), m_1_la)/(pow((tau_1_la+tau_2_la)/(tau_1_la+tau_1_la), m_1_la)+1)*1/(pow((tau_1_la+tau_2_la)/(tau_2_la+tau_2_la), m_2_la)+1);
              H_rv = pow((tau_1_rv+tau_2_rv)/(tau_1_rv+tau_1_rv), m_1_rv)/(pow((tau_1_rv+tau_2_rv)/(tau_1_rv+tau_1_rv), m_1_rv)+1)*1/(pow((tau_1_rv+tau_2_rv)/(tau_2_rv+tau_2_rv), m_2_rv)+1);
              H_ra = pow((tau_1_ra+tau_2_ra)/(tau_1_ra+tau_1_ra), m_1_ra)/(pow((tau_1_ra+tau_2_ra)/(tau_1_ra+tau_1_ra), m_1_ra)+1)*1/(pow((tau_1_ra+tau_2_ra)/(tau_2_ra+tau_2_ra), m_2_ra)+1);
              E_lv = (E_lv_max-E_lv_min)*g_1_lv/((g_1_lv+1)*(g_2_lv+1)*H_lv)+E_lv_min;
              g_1_lv = pow(mt/(tau_1_lv*T), m_1_lv);
              g_2_lv = pow(mt/(tau_2_lv*T), m_2_lv);
              E_la = (E_la_max-E_la_min)*g_1_la/((g_1_la+1)*(g_2_la+1)*H_la)+E_la_min;
              g_1_la = pow(mta/(tau_1_la*T), m_1_la);
              g_2_la = pow(mta/(tau_2_la*T), m_2_la);
              E_rv = (E_rv_max-E_rv_min)*g_1_rv/((g_1_rv+1)*(g_2_rv+1)*H_rv)+E_rv_min;
              g_1_rv = pow(mt/(tau_1_rv*T), m_1_rv);
              g_2_rv = pow(mt/(tau_2_rv*T), m_2_rv);
              E_ra = (E_ra_max-E_ra_min)*g_1_ra/((g_1_ra+1)*(g_2_ra+1)*H_ra)+E_ra_min;
              g_1_ra = pow(mta/(tau_1_ra*T), m_1_ra);
              g_2_ra = pow(mta/(tau_2_ra*T), m_2_ra);

              // Original model
        //       mt = t-T*floor(t/T);
        //       mta = t-t_onset_a*T-T*floor((t-t_onset_a*T)/T);

               //improvement for varying frequency
              der(mt) = 1;
              der(mta) = 1;

              when (time > pre(last_beat) + T) then
                reinit(mt, 0);
                last_beat = time;
              end when;

              when (time > last_beat + t_onset_a*T) then
                reinit(mta, 0);
              end when;



              //mta = t-t_onset_a*T-floor((t-t_onset_a*T));


              v_trv = noEvent(if u_ra >= u_rv and q_ra + q_ra_0 > 0 then
                      CQ_trv*sqrt(u_ra -u_rv)
                  else
                      0.0);
                       /*  u_ra < u_rv */

              v_puv = noEvent(if u_rv >= u_par and q_rv + q_rv_0 > 0 then
                      CQ_puv*sqrt(u_rv-u_par)
                  else
                      0.0);
                       /*  u_rv < u_par */

              v_miv = noEvent(if u_la >= u_lv and q_la + q_la_0 > 0 then
                      CQ_miv*sqrt(u_la-u_lv)
                  else
                      0.0);
                       /*  u_la < u_lv */

              v_aov = noEvent(if u_lv >= u_sas and q_lv + q_lv_0 > 0 then
                      CQ_aov*sqrt(u_lv-u_sas)
                  else
                      0.0);
                       /*  u_lv < u_sas */

              u_ra = E_ra*(q_ra-q_ra_0) + boundedTP_ra;
              u_rv = E_rv*(q_rv-q_rv_0) + boundedTP_rv;
              u_la = E_la*(q_la-q_la_0) + boundedTP_la;
              u_lv = E_lv*(q_lv-q_lv_0) + boundedTP_lv;


              der(q_ra) = v_sup_venacava+v_inf_venacava-v_trv;
              der(q_rv) = v_trv-v_puv;
              der(q_la) = v_pvn-v_miv;
              der(q_lv) = v_miv-v_aov;

        end HeartInputs;

        function pow
          input Real i;
          input Real j;
          output Real pou;
        algorithm
          pou :=i^j;
        end pow;

        model HeartComponent
          extends Auxiliary.HeartBase;
          HeartInputs Heart1(
            v_sup_venacava=sv.q,
            v_inf_venacava=0,
            u_sas=sa.pressure,
            u_par=pa.pressure,
            v_pvn=pv.q,
            t=time,
            T=1/HR_.y,
            thoracic_pressure = TP.y)
            annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

          Modelica.Blocks.Math.Gain TP(k=1) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=90,
                origin={0,-50})));
          Modelica.Blocks.Math.Gain HR_(k=1) annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-50,0})));
        equation

          sv.pressure = Heart1.u_ra;
          sa.q = -Heart1.v_aov;
          pv.pressure = Heart1.u_la;
          pa.q = -Heart1.v_puv;

          connect(TP.u, P0.y) annotation (Line(points={{-8.88178e-16,-62},{0,-62},{0,-75}},
                color={0,0,127}));
          connect(TP.u, thoracic_pressure_input) annotation (Line(points={{0,-62},{0,-82},
                  {0,-100},{-8,-100}}, color={0,0,127}));
          connect(HR0.y, HR_.u)
            annotation (Line(points={{-75,0},{-62,0}}, color={0,0,127}));
          connect(frequency_input, HR_.u)
            annotation (Line(points={{-106,0},{-62,0}}, color={0,0,127}));
          annotation (Icon(graphics={Text(
                  extent={{-100,60},{100,100}},
                  lineColor={0,0,0},
                  textString="ADAN_7af")}));
        end HeartComponent;
      end _7af7a4;

      package SystemicTissueParameters

      model SystemicTissueParameters
        import Physiolibrary.Types.*;
        parameter Pressure tissue_pressure=2666.4;
        parameter HydraulicResistance Ra_celiac_trunk_C116 = 8.65E+08 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_celiac_trunk_C116 = 1.62E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_celiac_trunk_C116 = 5.21E+05 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_celiac_trunk_C116 = 3.69E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_celiac_trunk_C116 = 2.85E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_renal_L166 = 9.58E+08 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_renal_L166 = 1.80E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_renal_L166 = 9.83E+05 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_renal_L166 = 3.33E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_renal_L166 = 2.57E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_renal_R178 = 9.61E+08 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_renal_R178 = 1.80E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_renal_R178 = 1.43E+06 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_renal_R178 = 3.32E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_renal_R178 = 2.56E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_internal_iliac_T1_R218 = 1.75E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_internal_iliac_T1_R218 = 3.28E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_internal_iliac_T1_R218 = 3.05E+06 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_internal_iliac_T1_R218 = 1.82E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_internal_iliac_T1_R218 = 1.41E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_profundus_T2_R224 = 1.53E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_profundus_T2_R224 = 2.87E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_profundus_T2_R224 = 1.73E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_profundus_T2_R224 = 2.08E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_profundus_T2_R224 = 1.61E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_anterior_tibial_T3_R230 = 1.20E+10 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_anterior_tibial_T3_R230 = 2.25E+09 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_anterior_tibial_T3_R230 = 9.49E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_anterior_tibial_T3_R230 = 2.65E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_anterior_tibial_T3_R230 = 2.05E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_posterior_tibial_T4_R236 = 1.03E+10 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_posterior_tibial_T4_R236 = 1.93E+09 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_posterior_tibial_T4_R236 = 8.47E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_posterior_tibial_T4_R236 = 3.10E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_posterior_tibial_T4_R236 = 2.39E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_internal_iliac_T1_L196 = 1.75E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_internal_iliac_T1_L196 = 3.29E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_internal_iliac_T1_L196 = 3.05E+06 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_internal_iliac_T1_L196 = 1.82E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_internal_iliac_T1_L196 = 1.40E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_profundus_T2_L202 = 1.53E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_profundus_T2_L202 = 2.87E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_profundus_T2_L202 = 1.73E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_profundus_T2_L202 = 2.08E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_profundus_T2_L202 = 1.61E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_anterior_tibial_T3_L208 = 1.20E+10 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_anterior_tibial_T3_L208 = 2.25E+09 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_anterior_tibial_T3_L208 = 9.49E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_anterior_tibial_T3_L208 = 2.65E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_anterior_tibial_T3_L208 = 2.05E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_posterior_tibial_T4_L214 = 1.03E+10 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_posterior_tibial_T4_L214 = 1.93E+09 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_posterior_tibial_T4_L214 = 8.47E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_posterior_tibial_T4_L214 = 3.10E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_posterior_tibial_T4_L214 = 2.39E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_ulnar_T2_R42 = 5.13E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_ulnar_T2_R42 = 9.61E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_ulnar_T2_R42 = 4.03E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_ulnar_T2_R42 = 6.22E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_ulnar_T2_R42 = 4.81E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_radial_T1_R44 = 6.11E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_radial_T1_R44 = 1.15E+09 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_radial_T1_R44 = 5.32E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_radial_T1_R44 = 5.21E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_radial_T1_R44 = 4.03E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_ulnar_T2_L90 = 5.21E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_ulnar_T2_L90 = 9.77E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_ulnar_T2_L90 = 4.03E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_ulnar_T2_L90 = 6.12E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_ulnar_T2_L90 = 4.73E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_radial_T1_L92 = 6.02E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_radial_T1_L92 = 1.13E+09 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_radial_T1_L92 = 5.32E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_radial_T1_L92 = 5.29E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_radial_T1_L92 = 4.09E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_internal_carotid_R8_C = 4.26E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_internal_carotid_R8_C = 7.99E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_internal_carotid_R8_C = 6.13E+06 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_internal_carotid_R8_C = 7.48E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_internal_carotid_R8_C = 5.78E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_external_carotid_T2_R26 = 4.19E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_external_carotid_T2_R26 = 7.86E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_external_carotid_T2_R26 = 3.97E+06 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_external_carotid_T2_R26 = 7.60E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_external_carotid_T2_R26 = 5.87E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_internal_carotid_L50_C = 4.28E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_internal_carotid_L50_C = 8.02E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_internal_carotid_L50_C = 6.13E+06 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_internal_carotid_L50_C = 7.45E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_internal_carotid_L50_C = 5.76E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_external_carotid_T2_L62 = 4.21E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_external_carotid_T2_L62 = 7.89E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_external_carotid_T2_L62 = 3.97E+06 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_external_carotid_T2_L62 = 7.57E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_external_carotid_T2_L62 = 5.85E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_vertebral_L2 = 4.70E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_vertebral_L2 = 8.80E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_vertebral_L2 = 3.93E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_vertebral_L2 = 6.79E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_vertebral_L2 = 5.25E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_vertebral_R272 = 4.68E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_vertebral_R272 = 8.77E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_vertebral_R272 = 3.94E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_vertebral_R272 = 6.81E-09 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_vertebral_R272 = 5.26E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_splachnic_tissue = 1.19E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_splachnic_tissue = 2.23E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_splachnic_tissue = 1.00E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_splachnic_tissue = 7.78E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_splachnic_tissue = 2.07E-04 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        parameter HydraulicResistance Ra_cardiac_tissue = 2.67E+09 annotation(Dialog(tab="Tissue parametrization", group="Arterioles resistance"));
        parameter HydraulicResistance Rv_cardiac_tissue = 5.01E+08 annotation(Dialog(tab="Tissue parametrization", group="Venules resistance"));
        parameter HydraulicInertance I_cardiac_tissue = 1.00E+07 annotation(Dialog(tab="Tissue parametrization", group="Inertances"));
        parameter HydraulicCompliance C_cardiac_tissue = 3.46E-08 annotation(Dialog(tab="Tissue parametrization", group="Compliances"));
        parameter Volume Zpv_cardiac_tissue = 9.22E-05 annotation(Dialog(tab="Tissue parametrization", group="Zero pressure volumes"));
        annotation(defaultComponentName = "tissueParameters");
      end SystemicTissueParameters;

        model SystemicTissueParameters_Calculated
          extends SystemicTissueParameters(
            Ra_celiac_trunk_C116 = (arterioles_pressure - tissue_pressure)/q_celiac_trunk_C116,
            Rv_celiac_trunk_C116 = (tissue_pressure - venules_pressure) /q_celiac_trunk_C116,
            C_celiac_trunk_C116 = stressed_volume/tissue_pressure*qf_celiac_trunk_C116,
            Zpv_celiac_trunk_C116 = total_zpv*qf_celiac_trunk_C116,
            Ra_renal_L166 = (arterioles_pressure - tissue_pressure)/q_renal_L166,
            Rv_renal_L166 = (tissue_pressure - venules_pressure) /q_renal_L166,
            C_renal_L166 = stressed_volume/tissue_pressure*qf_renal_L166,
            Zpv_renal_L166 = total_zpv*qf_renal_L166,
            Ra_renal_R178 = (arterioles_pressure - tissue_pressure)/q_renal_R178,
            Rv_renal_R178 = (tissue_pressure - venules_pressure) /q_renal_R178,
            C_renal_R178 = stressed_volume/tissue_pressure*qf_renal_R178,
            Zpv_renal_R178 = total_zpv*qf_renal_R178,
            Ra_internal_iliac_T1_R218 = (arterioles_pressure - tissue_pressure)/q_internal_iliac_T1_R218,
            Rv_internal_iliac_T1_R218 = (tissue_pressure - venules_pressure) /q_internal_iliac_T1_R218,
            C_internal_iliac_T1_R218 = stressed_volume/tissue_pressure*qf_internal_iliac_T1_R218,
            Zpv_internal_iliac_T1_R218 = total_zpv*qf_internal_iliac_T1_R218,
            Ra_profundus_T2_R224 = (arterioles_pressure - tissue_pressure)/q_profundus_T2_R224,
            Rv_profundus_T2_R224 = (tissue_pressure - venules_pressure) /q_profundus_T2_R224,
            C_profundus_T2_R224 = stressed_volume/tissue_pressure*qf_profundus_T2_R224,
            Zpv_profundus_T2_R224 = total_zpv*qf_profundus_T2_R224,
            Ra_anterior_tibial_T3_R230 = (arterioles_pressure - tissue_pressure)/q_anterior_tibial_T3_R230,
            Rv_anterior_tibial_T3_R230 = (tissue_pressure - venules_pressure) /q_anterior_tibial_T3_R230,
            C_anterior_tibial_T3_R230 = stressed_volume/tissue_pressure*qf_anterior_tibial_T3_R230,
            Zpv_anterior_tibial_T3_R230 = total_zpv*qf_anterior_tibial_T3_R230,
            Ra_posterior_tibial_T4_R236 = (arterioles_pressure - tissue_pressure)/q_posterior_tibial_T4_R236,
            Rv_posterior_tibial_T4_R236 = (tissue_pressure - venules_pressure) /q_posterior_tibial_T4_R236,
            C_posterior_tibial_T4_R236 = stressed_volume/tissue_pressure*qf_posterior_tibial_T4_R236,
            Zpv_posterior_tibial_T4_R236 = total_zpv*qf_posterior_tibial_T4_R236,
            Ra_internal_iliac_T1_L196 = (arterioles_pressure - tissue_pressure)/q_internal_iliac_T1_L196,
            Rv_internal_iliac_T1_L196 = (tissue_pressure - venules_pressure) /q_internal_iliac_T1_L196,
            C_internal_iliac_T1_L196 = stressed_volume/tissue_pressure*qf_internal_iliac_T1_L196,
            Zpv_internal_iliac_T1_L196 = total_zpv*qf_internal_iliac_T1_L196,
            Ra_profundus_T2_L202 = (arterioles_pressure - tissue_pressure)/q_profundus_T2_L202,
            Rv_profundus_T2_L202 = (tissue_pressure - venules_pressure) /q_profundus_T2_L202,
            C_profundus_T2_L202 = stressed_volume/tissue_pressure*qf_profundus_T2_L202,
            Zpv_profundus_T2_L202 = total_zpv*qf_profundus_T2_L202,
            Ra_anterior_tibial_T3_L208 = (arterioles_pressure - tissue_pressure)/q_anterior_tibial_T3_L208,
            Rv_anterior_tibial_T3_L208 = (tissue_pressure - venules_pressure) /q_anterior_tibial_T3_L208,
            C_anterior_tibial_T3_L208 = stressed_volume/tissue_pressure*qf_anterior_tibial_T3_L208,
            Zpv_anterior_tibial_T3_L208 = total_zpv*qf_anterior_tibial_T3_L208,
            Ra_posterior_tibial_T4_L214 = (arterioles_pressure - tissue_pressure)/q_posterior_tibial_T4_L214,
            Rv_posterior_tibial_T4_L214 = (tissue_pressure - venules_pressure) /q_posterior_tibial_T4_L214,
            C_posterior_tibial_T4_L214 = stressed_volume/tissue_pressure*qf_posterior_tibial_T4_L214,
            Zpv_posterior_tibial_T4_L214 = total_zpv*qf_posterior_tibial_T4_L214,
            Ra_ulnar_T2_R42 = (arterioles_pressure - tissue_pressure)/q_ulnar_T2_R42,
            Rv_ulnar_T2_R42 = (tissue_pressure - venules_pressure) /q_ulnar_T2_R42,
            C_ulnar_T2_R42 = stressed_volume/tissue_pressure*qf_ulnar_T2_R42,
            Zpv_ulnar_T2_R42 = total_zpv*qf_ulnar_T2_R42,
            Ra_radial_T1_R44 = (arterioles_pressure - tissue_pressure)/q_radial_T1_R44,
            Rv_radial_T1_R44 = (tissue_pressure - venules_pressure) /q_radial_T1_R44,
            C_radial_T1_R44 = stressed_volume/tissue_pressure*qf_radial_T1_R44,
            Zpv_radial_T1_R44 = total_zpv*qf_radial_T1_R44,
            Ra_ulnar_T2_L90 = (arterioles_pressure - tissue_pressure)/q_ulnar_T2_L90,
            Rv_ulnar_T2_L90 = (tissue_pressure - venules_pressure) /q_ulnar_T2_L90,
            C_ulnar_T2_L90 = stressed_volume/tissue_pressure*qf_ulnar_T2_L90,
            Zpv_ulnar_T2_L90 = total_zpv*qf_ulnar_T2_L90,
            Ra_radial_T1_L92 = (arterioles_pressure - tissue_pressure)/q_radial_T1_L92,
            Rv_radial_T1_L92 = (tissue_pressure - venules_pressure) /q_radial_T1_L92,
            C_radial_T1_L92 = stressed_volume/tissue_pressure*qf_radial_T1_L92,
            Zpv_radial_T1_L92 = total_zpv*qf_radial_T1_L92,
            Ra_internal_carotid_R8_C = (arterioles_pressure - tissue_pressure)/q_internal_carotid_R8_C,
            Rv_internal_carotid_R8_C = (tissue_pressure - venules_pressure) /q_internal_carotid_R8_C,
            C_internal_carotid_R8_C = stressed_volume/tissue_pressure*qf_internal_carotid_R8_C,
            Zpv_internal_carotid_R8_C = total_zpv*qf_internal_carotid_R8_C,
            Ra_external_carotid_T2_R26 = (arterioles_pressure - tissue_pressure)/q_external_carotid_T2_R26,
            Rv_external_carotid_T2_R26 = (tissue_pressure - venules_pressure) /q_external_carotid_T2_R26,
            C_external_carotid_T2_R26 = stressed_volume/tissue_pressure*qf_external_carotid_T2_R26,
            Zpv_external_carotid_T2_R26 = total_zpv*qf_external_carotid_T2_R26,
            Ra_internal_carotid_L50_C = (arterioles_pressure - tissue_pressure)/q_internal_carotid_L50_C,
            Rv_internal_carotid_L50_C = (tissue_pressure - venules_pressure) /q_internal_carotid_L50_C,
            C_internal_carotid_L50_C = stressed_volume/tissue_pressure*qf_internal_carotid_L50_C,
            Zpv_internal_carotid_L50_C = total_zpv*qf_internal_carotid_L50_C,
            Ra_external_carotid_T2_L62 = (arterioles_pressure - tissue_pressure)/q_external_carotid_T2_L62,
            Rv_external_carotid_T2_L62 = (tissue_pressure - venules_pressure) /q_external_carotid_T2_L62,
            C_external_carotid_T2_L62 = stressed_volume/tissue_pressure*qf_external_carotid_T2_L62,
            Zpv_external_carotid_T2_L62 = total_zpv*qf_external_carotid_T2_L62,
            Ra_vertebral_L2 = (arterioles_pressure - tissue_pressure)/q_vertebral_L2,
            Rv_vertebral_L2 = (tissue_pressure - venules_pressure) /q_vertebral_L2,
            C_vertebral_L2 = stressed_volume/tissue_pressure*qf_vertebral_L2,
            Zpv_vertebral_L2 = total_zpv*qf_vertebral_L2,
            Ra_vertebral_R272 = (arterioles_pressure - tissue_pressure)/q_vertebral_R272,
            Rv_vertebral_R272 = (tissue_pressure - venules_pressure) /q_vertebral_R272,
            C_vertebral_R272 = stressed_volume/tissue_pressure*qf_vertebral_R272,
            Zpv_vertebral_R272 = total_zpv*qf_vertebral_R272,
            Ra_splachnic_tissue = (arterioles_pressure - tissue_pressure)/q_splachnic_tissue,
            Rv_splachnic_tissue = (tissue_pressure - venules_pressure) /q_splachnic_tissue,
            C_splachnic_tissue = stressed_volume/tissue_pressure*qf_splachnic_tissue,
            Zpv_splachnic_tissue = total_zpv*qf_splachnic_tissue,
            Ra_cardiac_tissue = (arterioles_pressure - tissue_pressure)/q_cardiac_tissue,
            Rv_cardiac_tissue = (tissue_pressure - venules_pressure) /q_cardiac_tissue,
            C_cardiac_tissue = stressed_volume/tissue_pressure*qf_cardiac_tissue,
            Zpv_cardiac_tissue = total_zpv*qf_cardiac_tissue);
          import Physiolibrary.Types.*;

        parameter Pressure arterioles_pressure=13332;
        parameter Pressure venules_pressure=666.6;
        parameter Volume total_zpv=0.002304;
        parameter Volume stressed_volume=0.000795;
        parameter VolumeFlowRate cardiac_output=9.98e-5;

        parameter Fraction qf_celiac_trunk_C116=0.124          "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_renal_L166=0.112                 "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_renal_R178=0.111                 "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_internal_iliac_T1_R218=0.061     "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_profundus_T2_R224=0.07           "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_anterior_tibial_T3_R230=0.009    "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_posterior_tibial_T4_R236=0.01    "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_internal_iliac_T1_L196=0.061     "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_profundus_T2_L202=0.07           "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_anterior_tibial_T3_L208=0.009    "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_posterior_tibial_T4_L214=0.01    "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_ulnar_T2_R42=0.021               "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_radial_T1_R44=0.017              "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_ulnar_T2_L90=0.021               "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_radial_T1_L92=0.018              "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_internal_carotid_R8_C=0.025      "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_external_carotid_T2_R26=0.025    "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_internal_carotid_L50_C=0.025     "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_external_carotid_T2_L62=0.025    "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_vertebral_L2=0.023               "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_vertebral_R272=0.023             "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_splachnic_tissue=0.09            "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));
        parameter Fraction qf_cardiac_tissue=0.04              "Fraction of tissue size, estimated from flow fraction" annotation(Dialog(tab="General", group="Flow fractions"));

        parameter VolumeFlowRate q_celiac_trunk_C116 =        cardiac_output*qf_celiac_trunk_C116 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_renal_L166 =               cardiac_output*qf_renal_L166 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_renal_R178 =               cardiac_output*qf_renal_R178 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_internal_iliac_T1_R218 =   cardiac_output*qf_internal_iliac_T1_R218 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_profundus_T2_R224 =        cardiac_output*qf_profundus_T2_R224 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_anterior_tibial_T3_R230 =  cardiac_output*qf_anterior_tibial_T3_R230 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_posterior_tibial_T4_R236 = cardiac_output*qf_posterior_tibial_T4_R236 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_internal_iliac_T1_L196 =   cardiac_output*qf_internal_iliac_T1_L196 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_profundus_T2_L202 =        cardiac_output*qf_profundus_T2_L202 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_anterior_tibial_T3_L208 =  cardiac_output*qf_anterior_tibial_T3_L208 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_posterior_tibial_T4_L214 = cardiac_output*qf_posterior_tibial_T4_L214 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_ulnar_T2_R42 =             cardiac_output*qf_ulnar_T2_R42 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_radial_T1_R44 =            cardiac_output*qf_radial_T1_R44 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_ulnar_T2_L90 =             cardiac_output*qf_ulnar_T2_L90 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_radial_T1_L92 =            cardiac_output*qf_radial_T1_L92 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_internal_carotid_R8_C =    cardiac_output*qf_internal_carotid_R8_C "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_external_carotid_T2_R26 =  cardiac_output*qf_external_carotid_T2_R26 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_internal_carotid_L50_C =   cardiac_output*qf_internal_carotid_L50_C "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_external_carotid_T2_L62 =  cardiac_output*qf_external_carotid_T2_L62 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_vertebral_L2 =             cardiac_output*qf_vertebral_L2 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_vertebral_R272 =           cardiac_output*qf_vertebral_R272 "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_splachnic_tissue =         cardiac_output*qf_splachnic_tissue "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));
        parameter VolumeFlowRate q_cardiac_tissue =           cardiac_output*qf_cardiac_tissue "Tissue flow" annotation(Dialog(tab="Calculated", group="Flows"));

        end SystemicTissueParameters_Calculated;
         annotation (Documentation(info="<html>
<p>Generated tissue parameters</p>
</html>",       revisions="<html>
<p>This revision was generated at 2019-07-15 11:36:55.332316</p>
</html>"));
      end SystemicTissueParameters;

      partial model Systemic_interfaces
      //   parameter Real alphaC = 2.5;
      //   parameter Real alphaZPV = 2.5;
      //   inner Physiolibrary.Types.Fraction ZPV_effect = alphaZPV*(phi-0.25);
      //   inner Physiolibrary.Types.Fraction C_effect = alphaC*(phi-0.25);
      //   inner parameter Physiolibrary.Types.Fraction venous_diameter_correction = 1.5;
      //   inner parameter Physiolibrary.Types.Fraction C_fact = 1;
      //   inner parameter Physiolibrary.Types.Fraction cfactor = 1;
      //   inner parameter Real gamma = 1/2;
      //   inner parameter Real alpha = 5;// = T_pass/(T_pass + T_act_max);
      //   inner parameter Physiolibrary.Types.Fraction Ra_factor = 1 "Exponential factor affecting arterioles resistance";
      //
      //   inner parameter Physiolibrary.Types.Fraction phi0 = 0.25 "default value of phi. Also used for normalization";


        inner Modelica.SIunits.Angle Tilt;
        inner Physiolibrary.Types.Pressure thoracic_pressure;
        inner Physiolibrary.Types.Fraction phi "a systemic acitvation fraction, 1 being maximal possible. Normal resting is believed to be 1/4 of the maximum (0.25)";
        inner Physiolibrary.Types.Fraction Exercise;

        parameter Boolean UseThoracic_PressureInput = false annotation(choices(checkBox=true));
        parameter Boolean UsePhi_Input = false annotation(choices(checkBox=true));
        parameter Boolean UseTiltInput = false annotation(choices(checkBox=true));
        parameter Boolean UseExerciseInput = false annotation(choices(checkBox=true));


        replaceable model Systemic_artery_thoracic =
            ADAN_main.Components.Vessel_modules.pv_type_thoracic
                                                      constrainedby
          ADAN_main.Components.Vessel_modules.Interfaces.bg_vessel
                                                            annotation (choices(
         choice=ADAN_main.Components.Vessel_modules.pv_type_thoracic
          "No position calculations",
         choice=ADAN_main.Components.Vessel_modules.pv_type_thoracic_leveled
          "Position calculation"));
        replaceable model Systemic_artery =
            ADAN_main.Components.Vessel_modules.pv_type
          constrainedby
          ADAN_main.Components.Vessel_modules.Interfaces.bg_vessel
          annotation (choices(choice=ADAN_main.Components.Vessel_modules.pv_type
              "No position calculations", choice=ADAN_main.Components.Vessel_modules.pv_type_leveled
              "Position calculation"));
        replaceable model Systemic_vein =
            ADAN_main.Components.Vessel_modules.vp_type_tension_based constrainedby
          ADAN_main.Components.Vessel_modules.vp_type;

        replaceable model Systemic_tissue =
              ADAN_main.Components.Vessel_modules.systemic_tissue
          constrainedby
          ADAN_main.Components.Vessel_modules.Interfaces.systemic_tissue_base
                                                                 annotation (choices(
              choice=ADAN_main.Components.Vessel_modules.systemic_tissue
              "No position calculations", choice=ADAN_main.Components.Vessel_modules.systemic_tissue_leveled
              "Position calculation"));

        Modelica.Blocks.Interfaces.RealInput tilt_input = Tilt if UseTiltInput annotation (Placement(
            transformation(extent={{-320,-50},{-280,-10}}),
                                                          iconTransformation(extent={{26,-100},
                  {66,-60}})));
        Physiolibrary.Types.RealIO.PressureInput thoracic_pressure_input = thoracic_pressure if UseThoracic_PressureInput annotation (Placement(
            transformation(extent={{-320,-90},{-280,-50}}),
                                                          iconTransformation(extent={{-40,
                  -100},{0,-60}})));
        Physiolibrary.Types.RealIO.FractionInput phi_input = phi if UsePhi_Input annotation (Placement(
            transformation(extent={{-280,-50},{-240,-10}}),
                                                          iconTransformation(extent={{100,
                  -100},{140,-60}})));

      //   Physiolibrary.Types.RealIO.PressureOutput  p = internal_carotid_R8_A.u_in
      //                                                 annotation (Placement(
      //         transformation(extent={{-42,-92},{-22,-72}}),
      //                                                     iconTransformation(extent={{130,162},
      //             {150,182}})));


        Physiolibrary.Types.RealIO.PressureInput exercise_input=Exercise if            UseExerciseInput
          annotation (Placement(transformation(extent={{-320,-90},{-280,-50}}),
              iconTransformation(extent={{-40,60},{0,100}})));
        outer Settings settings annotation (Placement(transformation(extent={{
                  -318,180},{-298,200}})));

      equation
        if not UseThoracic_PressureInput then
          thoracic_pressure = 0;
        end if;
        if not UsePhi_Input then
          phi = settings.phi0;
        end if;
        if not UseTiltInput then
          Tilt = 0;
        end if;
        if not UseExerciseInput then
          Exercise = 0;
        end if;
          annotation (Placement(transformation(extent={{-100,80},{-80,100}})),
                         choices(choice=ADAN_main.Components.Vessel_modules.vp_type
                                                                                   "No position calculations",
             choice=ADAN_main.Components.Vessel_modules.vp_type_leveled "Position calculation",
             choice=ADAN_main.Components.Vessel_modules.vp_type_phi_sensitive "Phi_sensitive"),
                    Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Systemic_interfaces;

      partial model Systemic_base
        extends Systemic_interfaces;

        ADAN_main.Components.AdanVenousRed._b580e.Parameters_Venous_cellml.Parameters_Systemic
          Parameters_Systemic1
          annotation (Placement(transformation(extent={{-96,-87},{-76,-82}})));
        replaceable ADAN_main.Components.Vessel_modules.vv_type_thoracic
          ascending_aorta_A constrainedby
          ADAN_main.Components.Vessel_modules.Interfaces.bg_vessel(
          l=Parameters_Systemic1.l_ascending_aorta_A,
          E=Parameters_Systemic1.E_ascending_aorta_A,
          r=Parameters_Systemic1.r_ascending_aorta_A) annotation (Placement(
              transformation(extent={{-263,85},{-243,90}})),
            __Dymola_choicesAllMatching=true);

        Systemic_artery_thoracic ascending_aorta_B(
            l = Parameters_Systemic1.l_ascending_aorta_B,
            E = Parameters_Systemic1.E_ascending_aorta_B,
            r = Parameters_Systemic1.r_ascending_aorta_B)
        annotation (Placement(transformation(extent={{-238,85},{-218,90}})));

        Systemic_artery_thoracic ascending_aorta_C(
            l = Parameters_Systemic1.l_ascending_aorta_C,
            E = Parameters_Systemic1.E_ascending_aorta_C,
            r = Parameters_Systemic1.r_ascending_aorta_C)
        annotation (Placement(transformation(extent={{-213,85},{-193,90}})));

        Systemic_artery_thoracic ascending_aorta_D(
            l = Parameters_Systemic1.l_ascending_aorta_D,
            E = Parameters_Systemic1.E_ascending_aorta_D,
            r = Parameters_Systemic1.r_ascending_aorta_D)
        annotation (Placement(transformation(extent={{-188,85},{-168,90}})));

          Systemic_artery_thoracic aortic_arch_C2(
            l=Parameters_Systemic1.l_aortic_arch_C2,
            E=Parameters_Systemic1.E_aortic_arch_C2,
            r=Parameters_Systemic1.r_aortic_arch_C2)
            annotation (Placement(transformation(extent={{-163,85},{-143,90}})));

          Systemic_artery_thoracic
            brachiocephalic_trunk_C4(
            l=Parameters_Systemic1.l_brachiocephalic_trunk_C4,
            E=Parameters_Systemic1.E_brachiocephalic_trunk_C4,
            r=Parameters_Systemic1.r_brachiocephalic_trunk_C4)
            annotation (Placement(transformation(extent={{-132,163},{-112,168}})));

          Systemic_artery_thoracic
            aortic_arch_C46(
            l=Parameters_Systemic1.l_aortic_arch_C46,
            E=Parameters_Systemic1.E_aortic_arch_C46,
            r=Parameters_Systemic1.r_aortic_arch_C46)
            annotation (Placement(transformation(extent={{-137,85},{-117,90}})));

          Systemic_artery_thoracic aortic_arch_C64(
            l=Parameters_Systemic1.l_aortic_arch_C64,
            E=Parameters_Systemic1.E_aortic_arch_C64,
            r=Parameters_Systemic1.r_aortic_arch_C64)
            annotation (Placement(transformation(extent={{-112,85},{-92,90}})));

        Systemic_artery_thoracic aortic_arch_C94(
            l = Parameters_Systemic1.l_aortic_arch_C94,
            E = Parameters_Systemic1.E_aortic_arch_C94,
            r = Parameters_Systemic1.r_aortic_arch_C94)
        annotation (Placement(transformation(extent={{-91,65},{-111,70}})));

        Systemic_artery_thoracic thoracic_aorta_C96(
            l = Parameters_Systemic1.l_thoracic_aorta_C96,
            E = Parameters_Systemic1.E_thoracic_aorta_C96,
            r = Parameters_Systemic1.r_thoracic_aorta_C96)
        annotation (Placement(transformation(extent={{-118,65},{-138,70}})));

        Systemic_artery_thoracic thoracic_aorta_C100(
            l = Parameters_Systemic1.l_thoracic_aorta_C100,
            E = Parameters_Systemic1.E_thoracic_aorta_C100,
            r = Parameters_Systemic1.r_thoracic_aorta_C100)
        annotation (Placement(transformation(extent={{-143,65},{-163,70}})));

        Systemic_artery_thoracic thoracic_aorta_C104(
            l = Parameters_Systemic1.l_thoracic_aorta_C104,
            E = Parameters_Systemic1.E_thoracic_aorta_C104,
            r = Parameters_Systemic1.r_thoracic_aorta_C104)
        annotation (Placement(transformation(extent={{-168,65},{-188,70}})));

        Systemic_artery_thoracic thoracic_aorta_C108(
            l = Parameters_Systemic1.l_thoracic_aorta_C108,
            E = Parameters_Systemic1.E_thoracic_aorta_C108,
            r = Parameters_Systemic1.r_thoracic_aorta_C108)
        annotation (Placement(transformation(extent={{-193,65},{-213,70}})));

        Systemic_artery_thoracic thoracic_aorta_C112(
            l = Parameters_Systemic1.l_thoracic_aorta_C112,
            E = Parameters_Systemic1.E_thoracic_aorta_C112,
            r = Parameters_Systemic1.r_thoracic_aorta_C112)
        annotation (Placement(transformation(extent={{-218,65},{-238,70}})));

          Systemic_artery abdominal_aorta_C114(
            l=Parameters_Systemic1.l_abdominal_aorta_C114,
            E=Parameters_Systemic1.E_abdominal_aorta_C114,
            r=Parameters_Systemic1.r_abdominal_aorta_C114)
            annotation (Placement(transformation(extent={{-307,-5},{-287,0}})));
        Systemic_artery abdominal_aorta_C136(
            l = Parameters_Systemic1.l_abdominal_aorta_C136,
            E = Parameters_Systemic1.E_abdominal_aorta_C136,
            r = Parameters_Systemic1.r_abdominal_aorta_C136)
        annotation (Placement(transformation(extent={{-280,-5},{-260,0}})));
          Systemic_artery abdominal_aorta_C164(
            l=Parameters_Systemic1.l_abdominal_aorta_C164,
            E=Parameters_Systemic1.E_abdominal_aorta_C164,
            r=Parameters_Systemic1.r_abdominal_aorta_C164)
            annotation (Placement(transformation(extent={{-255,-5},{-235,0}})));
          Systemic_artery abdominal_aorta_C176(
            l=Parameters_Systemic1.l_abdominal_aorta_C176,
            E=Parameters_Systemic1.E_abdominal_aorta_C176,
            r=Parameters_Systemic1.r_abdominal_aorta_C176)
            annotation (Placement(transformation(extent={{-230,-5},{-210,0}})));
        Systemic_artery abdominal_aorta_C188(
            l = Parameters_Systemic1.l_abdominal_aorta_C188,
            E = Parameters_Systemic1.E_abdominal_aorta_C188,
            r = Parameters_Systemic1.r_abdominal_aorta_C188)
        annotation (Placement(transformation(extent={{-203,-5},{-183,0}})));
          Systemic_artery abdominal_aorta_C192(
            l=Parameters_Systemic1.l_abdominal_aorta_C192,
            E=Parameters_Systemic1.E_abdominal_aorta_C192,
            r=Parameters_Systemic1.r_abdominal_aorta_C192)
            annotation (Placement(transformation(extent={{-178,-5},{-158,0}})));
        Systemic_tissue celiac_trunk_C116(
            Ra = tissueParameters.Ra_celiac_trunk_C116,
            Rv = tissueParameters.Rv_celiac_trunk_C116,
            I = tissueParameters.I_celiac_trunk_C116,
            C = tissueParameters.C_celiac_trunk_C116,
            zpv =  tissueParameters.Zpv_celiac_trunk_C116,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,57},{55,62}})));
        Systemic_tissue renal_L166(
            Ra = tissueParameters.Ra_renal_L166,
            Rv = tissueParameters.Rv_renal_L166,
            I = tissueParameters.I_renal_L166,
            C = tissueParameters.C_renal_L166,
            zpv =  tissueParameters.Zpv_renal_L166,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,-67},{55,-62}})));
        Systemic_tissue renal_R178(
            Ra = tissueParameters.Ra_renal_R178,
            Rv = tissueParameters.Rv_renal_R178,
            I = tissueParameters.I_renal_R178,
            C = tissueParameters.C_renal_R178,
            zpv =  tissueParameters.Zpv_renal_R178,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,37},{55,42}})));
          Systemic_artery common_iliac_R216(
            l=Parameters_Systemic1.l_common_iliac_R216,
            E=Parameters_Systemic1.E_common_iliac_R216,
            r=Parameters_Systemic1.r_common_iliac_R216)
            annotation (Placement(transformation(extent={{-148,-5},{-128,0}})));
        Systemic_tissue internal_iliac_T1_R218(
            Ra = tissueParameters.Ra_internal_iliac_T1_R218,
            Rv = tissueParameters.Rv_internal_iliac_T1_R218,
            I = tissueParameters.I_internal_iliac_T1_R218,
            C = tissueParameters.C_internal_iliac_T1_R218,
            zpv =  tissueParameters.Zpv_internal_iliac_T1_R218,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,27},{55,32}})));
        Systemic_artery external_iliac_R220(
            l = Parameters_Systemic1.l_external_iliac_R220,
            E = Parameters_Systemic1.E_external_iliac_R220,
            r = Parameters_Systemic1.r_external_iliac_R220)
        annotation (Placement(transformation(extent={{-120,-5},{-100,0}})));
          Systemic_artery femoral_R222(
            l=Parameters_Systemic1.l_femoral_R222,
            E=Parameters_Systemic1.E_femoral_R222,
            r=Parameters_Systemic1.r_femoral_R222)
            annotation (Placement(transformation(extent={{-95,-5},{-75,0}})));
        Systemic_tissue profundus_T2_R224(
            Ra = tissueParameters.Ra_profundus_T2_R224,
            Rv = tissueParameters.Rv_profundus_T2_R224,
            I = tissueParameters.I_profundus_T2_R224,
            C = tissueParameters.C_profundus_T2_R224,
            zpv =  tissueParameters.Zpv_profundus_T2_R224,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,17},{55,22}})));
        Systemic_artery femoral_R226(
            l = Parameters_Systemic1.l_femoral_R226,
            E = Parameters_Systemic1.E_femoral_R226,
            r = Parameters_Systemic1.r_femoral_R226)
        annotation (Placement(transformation(extent={{-67,-5},{-47,0}})));
          Systemic_artery popliteal_R228(
            l=Parameters_Systemic1.l_popliteal_R228,
            E=Parameters_Systemic1.E_popliteal_R228,
            r=Parameters_Systemic1.r_popliteal_R228)
            annotation (Placement(transformation(extent={{-42,-5},{-22,0}})));
        Systemic_tissue anterior_tibial_T3_R230(
            Ra = tissueParameters.Ra_anterior_tibial_T3_R230,
            Rv = tissueParameters.Rv_anterior_tibial_T3_R230,
            I = tissueParameters.I_anterior_tibial_T3_R230,
            C = tissueParameters.C_anterior_tibial_T3_R230,
            zpv =  tissueParameters.Zpv_anterior_tibial_T3_R230,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,7},{55,12}})));
        Systemic_artery popliteal_R232(
            l = Parameters_Systemic1.l_popliteal_R232,
            E = Parameters_Systemic1.E_popliteal_R232,
            r = Parameters_Systemic1.r_popliteal_R232)
        annotation (Placement(transformation(extent={{-16,-5},{4,0}})));
        Systemic_artery tibiofibular_trunk_R234(
            l = Parameters_Systemic1.l_tibiofibular_trunk_R234,
            E = Parameters_Systemic1.E_tibiofibular_trunk_R234,
            r = Parameters_Systemic1.r_tibiofibular_trunk_R234)
        annotation (Placement(transformation(extent={{9,-5},{29,0}})));
        Systemic_tissue posterior_tibial_T4_R236(
            Ra = tissueParameters.Ra_posterior_tibial_T4_R236,
            Rv = tissueParameters.Rv_posterior_tibial_T4_R236,
            I = tissueParameters.I_posterior_tibial_T4_R236,
            C = tissueParameters.C_posterior_tibial_T4_R236,
            zpv =  tissueParameters.Zpv_posterior_tibial_T4_R236,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,-5},{55,0}})));
          Systemic_artery common_iliac_L194(
            l=Parameters_Systemic1.l_common_iliac_L194,
            E=Parameters_Systemic1.E_common_iliac_L194,
            r=Parameters_Systemic1.r_common_iliac_L194)
            annotation (Placement(transformation(extent={{-147,-23},{-127,-18}})));
        Systemic_tissue internal_iliac_T1_L196(
            Ra = tissueParameters.Ra_internal_iliac_T1_L196,
            Rv = tissueParameters.Rv_internal_iliac_T1_L196,
            I = tissueParameters.I_internal_iliac_T1_L196,
            C = tissueParameters.C_internal_iliac_T1_L196,
            zpv =  tissueParameters.Zpv_internal_iliac_T1_L196,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,-55},{55,-50}})));
        Systemic_artery external_iliac_L198(
            l = Parameters_Systemic1.l_external_iliac_L198,
            E = Parameters_Systemic1.E_external_iliac_L198,
            r = Parameters_Systemic1.r_external_iliac_L198)
        annotation (Placement(transformation(extent={{-121,-23},{-101,-18}})));
          Systemic_artery femoral_L200(
            l=Parameters_Systemic1.l_femoral_L200,
            E=Parameters_Systemic1.E_femoral_L200,
            r=Parameters_Systemic1.r_femoral_L200)
            annotation (Placement(transformation(extent={{-94,-23},{-74,-18}})));
        Systemic_tissue profundus_T2_L202(
            Ra = tissueParameters.Ra_profundus_T2_L202,
            Rv = tissueParameters.Rv_profundus_T2_L202,
            I = tissueParameters.I_profundus_T2_L202,
            C = tissueParameters.C_profundus_T2_L202,
            zpv =  tissueParameters.Zpv_profundus_T2_L202,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,-43},{55,-38}})));
        Systemic_artery femoral_L204(
            l = Parameters_Systemic1.l_femoral_L204,
            E = Parameters_Systemic1.E_femoral_L204,
            r = Parameters_Systemic1.r_femoral_L204)
        annotation (Placement(transformation(extent={{-68,-23},{-48,-18}})));
          Systemic_artery popliteal_L206(
            l=Parameters_Systemic1.l_popliteal_L206,
            E=Parameters_Systemic1.E_popliteal_L206,
            r=Parameters_Systemic1.r_popliteal_L206)
            annotation (Placement(transformation(extent={{-41,-23},{-21,-18}})));
        Systemic_tissue anterior_tibial_T3_L208(
            Ra = tissueParameters.Ra_anterior_tibial_T3_L208,
            Rv = tissueParameters.Rv_anterior_tibial_T3_L208,
            I = tissueParameters.I_anterior_tibial_T3_L208,
            C = tissueParameters.C_anterior_tibial_T3_L208,
            zpv =  tissueParameters.Zpv_anterior_tibial_T3_L208,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,-33},{55,-28}})));
        Systemic_artery popliteal_L210(
            l = Parameters_Systemic1.l_popliteal_L210,
            E = Parameters_Systemic1.E_popliteal_L210,
            r = Parameters_Systemic1.r_popliteal_L210)
        annotation (Placement(transformation(extent={{-16,-23},{4,-18}})));
        Systemic_artery tibiofibular_trunk_L212(
            l = Parameters_Systemic1.l_tibiofibular_trunk_L212,
            E = Parameters_Systemic1.E_tibiofibular_trunk_L212,
            r = Parameters_Systemic1.r_tibiofibular_trunk_L212)
        annotation (Placement(transformation(extent={{9,-23},{29,-18}})));
        Systemic_tissue posterior_tibial_T4_L214(
            Ra = tissueParameters.Ra_posterior_tibial_T4_L214,
            Rv = tissueParameters.Rv_posterior_tibial_T4_L214,
            I = tissueParameters.I_posterior_tibial_T4_L214,
            C = tissueParameters.C_posterior_tibial_T4_L214,
            zpv =  tissueParameters.Zpv_posterior_tibial_T4_L214,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,-23},{55,-18}})));
          Systemic_artery subclavian_R28(
            l=Parameters_Systemic1.l_subclavian_R28,
            E=Parameters_Systemic1.E_subclavian_R28,
            r=Parameters_Systemic1.r_subclavian_R28)
            annotation (Placement(transformation(extent={{-88,163},{-68,168}})));
        Systemic_artery subclavian_R30(
            l = Parameters_Systemic1.l_subclavian_R30,
            E = Parameters_Systemic1.E_subclavian_R30,
            r = Parameters_Systemic1.r_subclavian_R30)
        annotation (Placement(transformation(extent={{-63,163},{-43,168}})));
        Systemic_artery axillary_R32(
            l = Parameters_Systemic1.l_axillary_R32,
            E = Parameters_Systemic1.E_axillary_R32,
            r = Parameters_Systemic1.r_axillary_R32)
        annotation (Placement(transformation(extent={{-40,163},{-20,168}})));
          Systemic_artery brachial_R34(
            l=Parameters_Systemic1.l_brachial_R34,
            E=Parameters_Systemic1.E_brachial_R34,
            r=Parameters_Systemic1.r_brachial_R34)
            annotation (Placement(transformation(extent={{-17,163},{3,168}})));
        Systemic_artery ulnar_T2_R36(
            l = Parameters_Systemic1.l_ulnar_T2_R36,
            E = Parameters_Systemic1.E_ulnar_T2_R36,
            r = Parameters_Systemic1.r_ulnar_T2_R36)
        annotation (Placement(transformation(extent={{10,163},{30,168}})));
        Systemic_tissue ulnar_T2_R42(
            Ra = tissueParameters.Ra_ulnar_T2_R42,
            Rv = tissueParameters.Rv_ulnar_T2_R42,
            I = tissueParameters.I_ulnar_T2_R42,
            C = tissueParameters.C_ulnar_T2_R42,
            zpv =  tissueParameters.Zpv_ulnar_T2_R42,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,163},{55,168}})));
        Systemic_tissue radial_T1_R44(
            Ra = tissueParameters.Ra_radial_T1_R44,
            Rv = tissueParameters.Rv_radial_T1_R44,
            I = tissueParameters.I_radial_T1_R44,
            C = tissueParameters.C_radial_T1_R44,
            zpv =  tissueParameters.Zpv_radial_T1_R44,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,153},{55,158}})));
          Systemic_artery subclavian_L66(
            l=Parameters_Systemic1.l_subclavian_L66,
            E=Parameters_Systemic1.E_subclavian_L66,
            r=Parameters_Systemic1.r_subclavian_L66)
            annotation (Placement(transformation(extent={{-87,109},{-67,114}})));
        Systemic_artery subclavian_L78(
            l = Parameters_Systemic1.l_subclavian_L78,
            E = Parameters_Systemic1.E_subclavian_L78,
            r = Parameters_Systemic1.r_subclavian_L78)
        annotation (Placement(transformation(extent={{-62,109},{-42,114}})));
        Systemic_artery axillary_L80(
            l = Parameters_Systemic1.l_axillary_L80,
            E = Parameters_Systemic1.E_axillary_L80,
            r = Parameters_Systemic1.r_axillary_L80)
        annotation (Placement(transformation(extent={{-39,109},{-19,114}})));
          Systemic_artery brachial_L82(
            l=Parameters_Systemic1.l_brachial_L82,
            E=Parameters_Systemic1.E_brachial_L82,
            r=Parameters_Systemic1.r_brachial_L82)
            annotation (Placement(transformation(extent={{-16,109},{4,114}})));
        Systemic_artery ulnar_T2_L84(
            l = Parameters_Systemic1.l_ulnar_T2_L84,
            E = Parameters_Systemic1.E_ulnar_T2_L84,
            r = Parameters_Systemic1.r_ulnar_T2_L84)
        annotation (Placement(transformation(extent={{9,109},{29,114}})));
        Systemic_tissue ulnar_T2_L90(
            Ra = tissueParameters.Ra_ulnar_T2_L90,
            Rv = tissueParameters.Rv_ulnar_T2_L90,
            I = tissueParameters.I_ulnar_T2_L90,
            C = tissueParameters.C_ulnar_T2_L90,
            zpv =  tissueParameters.Zpv_ulnar_T2_L90,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,109},{55,114}})));
        Systemic_tissue radial_T1_L92(
            Ra = tissueParameters.Ra_radial_T1_L92,
            Rv = tissueParameters.Rv_radial_T1_L92,
            I = tissueParameters.I_radial_T1_L92,
            C = tissueParameters.C_radial_T1_L92,
            zpv =  tissueParameters.Zpv_radial_T1_L92,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,97},{55,102}})));
        Systemic_artery common_carotid_R6_A(
            l = Parameters_Systemic1.l_common_carotid_R6_A,
            E = Parameters_Systemic1.E_common_carotid_R6_A,
            r = Parameters_Systemic1.r_common_carotid_R6_A)
        annotation (Placement(transformation(extent={{-88,187},{-68,192}})));
        Systemic_artery common_carotid_R6_B(
            l = Parameters_Systemic1.l_common_carotid_R6_B,
            E = Parameters_Systemic1.E_common_carotid_R6_B,
            r = Parameters_Systemic1.r_common_carotid_R6_B)
        annotation (Placement(transformation(extent={{-63,187},{-43,192}})));
          Systemic_artery common_carotid_R6_C(
            l=Parameters_Systemic1.l_common_carotid_R6_C,
            E=Parameters_Systemic1.E_common_carotid_R6_C,
            r=Parameters_Systemic1.r_common_carotid_R6_C)
            annotation (Placement(transformation(extent={{-40,187},{-20,192}})));
        replaceable Systemic_artery internal_carotid_R8_A(
            l = Parameters_Systemic1.l_internal_carotid_R8_A,
            E = Parameters_Systemic1.E_internal_carotid_R8_A,
            r = Parameters_Systemic1.r_internal_carotid_R8_A)
        annotation (Placement(transformation(extent={{-17,187},{3,192}})));
        Systemic_artery internal_carotid_R8_B(
            l = Parameters_Systemic1.l_internal_carotid_R8_B,
            E = Parameters_Systemic1.E_internal_carotid_R8_B,
            r = Parameters_Systemic1.r_internal_carotid_R8_B)
        annotation (Placement(transformation(extent={{10,187},{30,192}})));
        Systemic_tissue internal_carotid_R8_C(
            Ra = tissueParameters.Ra_internal_carotid_R8_C,
            Rv = tissueParameters.Rv_internal_carotid_R8_C,
            I = tissueParameters.I_internal_carotid_R8_C,
            C = tissueParameters.C_internal_carotid_R8_C,
            zpv =  tissueParameters.Zpv_internal_carotid_R8_C,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,187},{55,192}})));
        Systemic_tissue external_carotid_T2_R26(
            Ra = tissueParameters.Ra_external_carotid_T2_R26,
            Rv = tissueParameters.Rv_external_carotid_T2_R26,
            I = tissueParameters.I_external_carotid_T2_R26,
            C = tissueParameters.C_external_carotid_T2_R26,
            zpv =  tissueParameters.Zpv_external_carotid_T2_R26,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,175},{55,180}})));
        Systemic_artery common_carotid_L48_A(
            l = Parameters_Systemic1.l_common_carotid_L48_A,
            E = Parameters_Systemic1.E_common_carotid_L48_A,
            r = Parameters_Systemic1.r_common_carotid_L48_A)
        annotation (Placement(transformation(extent={{-111,119},{-91,124}})));
        Systemic_artery common_carotid_L48_B(
            l = Parameters_Systemic1.l_common_carotid_L48_B,
            E = Parameters_Systemic1.E_common_carotid_L48_B,
            r = Parameters_Systemic1.r_common_carotid_L48_B)
        annotation (Placement(transformation(extent={{-86,119},{-66,124}})));
        Systemic_artery common_carotid_L48_C(
            l = Parameters_Systemic1.l_common_carotid_L48_C,
            E = Parameters_Systemic1.E_common_carotid_L48_C,
            r = Parameters_Systemic1.r_common_carotid_L48_C)
        annotation (Placement(transformation(extent={{-63,119},{-43,124}})));
          Systemic_artery common_carotid_L48_D(
            l=Parameters_Systemic1.l_common_carotid_L48_D,
            E=Parameters_Systemic1.E_common_carotid_L48_D,
            r=Parameters_Systemic1.r_common_carotid_L48_D)
            annotation (Placement(transformation(extent={{-38,119},{-18,124}})));
        Systemic_artery internal_carotid_L50_A(
            l = Parameters_Systemic1.l_internal_carotid_L50_A,
            E = Parameters_Systemic1.E_internal_carotid_L50_A,
            r = Parameters_Systemic1.r_internal_carotid_L50_A)
        annotation (Placement(transformation(extent={{-15,127},{5,132}})));
        Systemic_artery internal_carotid_L50_B(
            l = Parameters_Systemic1.l_internal_carotid_L50_B,
            E = Parameters_Systemic1.E_internal_carotid_L50_B,
            r = Parameters_Systemic1.r_internal_carotid_L50_B)
        annotation (Placement(transformation(extent={{10,127},{30,132}})));
        Systemic_tissue internal_carotid_L50_C(
            Ra = tissueParameters.Ra_internal_carotid_L50_C,
            Rv = tissueParameters.Rv_internal_carotid_L50_C,
            I = tissueParameters.I_internal_carotid_L50_C,
            C = tissueParameters.C_internal_carotid_L50_C,
            zpv =  tissueParameters.Zpv_internal_carotid_L50_C,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,127},{55,132}})));
        Systemic_tissue external_carotid_T2_L62(
            Ra = tissueParameters.Ra_external_carotid_T2_L62,
            Rv = tissueParameters.Rv_external_carotid_T2_L62,
            I = tissueParameters.I_external_carotid_T2_L62,
            C = tissueParameters.C_external_carotid_T2_L62,
            zpv =  tissueParameters.Zpv_external_carotid_T2_L62,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,119},{55,124}})));
        Systemic_tissue vertebral_L2(
            Ra = tissueParameters.Ra_vertebral_L2,
            Rv = tissueParameters.Rv_vertebral_L2,
            I = tissueParameters.I_vertebral_L2,
            C = tissueParameters.C_vertebral_L2,
            zpv =  tissueParameters.Zpv_vertebral_L2,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,87},{55,92}})));
        Systemic_tissue vertebral_R272(
            Ra = tissueParameters.Ra_vertebral_R272,
            Rv = tissueParameters.Rv_vertebral_R272,
            I = tissueParameters.I_vertebral_R272,
            C = tissueParameters.C_vertebral_R272,
            zpv =  tissueParameters.Zpv_vertebral_R272,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,143},{55,148}})));
        Systemic_vein superior_vena_cava_C2(
            l = Parameters_Venous1.l_superior_vena_cava_C2,
            E = Parameters_Venous1.E_superior_vena_cava_C2,
            r = Parameters_Venous1.r_superior_vena_cava_C2)
        annotation (Placement(transformation(extent={{297,163},{317,168}})));
          Systemic_vein superior_vena_cava_C88(
            l=Parameters_Venous1.l_superior_vena_cava_C88,
            E=Parameters_Venous1.E_superior_vena_cava_C88,
            r=Parameters_Venous1.r_superior_vena_cava_C88)
            annotation (Placement(transformation(extent={{272,163},{292,168}})));
          Systemic_vein inferior_vena_cava_C8(
            l=Parameters_Venous1.l_inferior_vena_cava_C8,
            E=Parameters_Venous1.E_inferior_vena_cava_C8,
            r=Parameters_Venous1.r_inferior_vena_cava_C8)
            annotation (Placement(transformation(extent={{409,-5},{429,0}})));
        Systemic_vein hepatic_vein_T1_C10(
            l = Parameters_Venous1.l_hepatic_vein_T1_C10,
            E = Parameters_Venous1.E_hepatic_vein_T1_C10,
            r = Parameters_Venous1.r_hepatic_vein_T1_C10)
        annotation (Placement(transformation(extent={{60,57},{80,62}})));
        Systemic_vein inferior_vena_cava_C12(
            l = Parameters_Venous1.l_inferior_vena_cava_C12,
            E = Parameters_Venous1.E_inferior_vena_cava_C12,
            r = Parameters_Venous1.r_inferior_vena_cava_C12)
        annotation (Placement(transformation(extent={{383,-5},{403,0}})));
          Systemic_vein inferior_vena_cava_C16(
            l=Parameters_Venous1.l_inferior_vena_cava_C16,
            E=Parameters_Venous1.E_inferior_vena_cava_C16,
            r=Parameters_Venous1.r_inferior_vena_cava_C16)
            annotation (Placement(transformation(extent={{360,-5},{380,0}})));
        Systemic_vein renal_vein_T1_R18(
            l = Parameters_Venous1.l_renal_vein_T1_R18,
            E = Parameters_Venous1.E_renal_vein_T1_R18,
            r = Parameters_Venous1.r_renal_vein_T1_R18)
        annotation (Placement(transformation(extent={{60,37},{80,42}})));
          Systemic_vein inferior_vena_cava_C20(
            l=Parameters_Venous1.l_inferior_vena_cava_C20,
            E=Parameters_Venous1.E_inferior_vena_cava_C20,
            r=Parameters_Venous1.r_inferior_vena_cava_C20)
            annotation (Placement(transformation(extent={{336,-5},{356,0}})));
        Systemic_vein renal_vein_T1_L22(
            l = Parameters_Venous1.l_renal_vein_T1_L22,
            E = Parameters_Venous1.E_renal_vein_T1_L22,
            r = Parameters_Venous1.r_renal_vein_T1_L22)
        annotation (Placement(transformation(extent={{60,-67},{80,-62}})));
          Systemic_vein inferior_vena_cava_C24(
            l=Parameters_Venous1.l_inferior_vena_cava_C24,
            E=Parameters_Venous1.E_inferior_vena_cava_C24,
            r=Parameters_Venous1.r_inferior_vena_cava_C24)
            annotation (Placement(transformation(extent={{312,-5},{332,0}})));
        Systemic_vein common_iliac_vein_L56(
            l = Parameters_Venous1.l_common_iliac_vein_L56,
            E = Parameters_Venous1.E_common_iliac_vein_L56,
            r = Parameters_Venous1.r_common_iliac_vein_L56)
        annotation (Placement(transformation(extent={{287,-23},{307,-18}})));
        Systemic_vein common_iliac_vein_R26(
            l = Parameters_Venous1.l_common_iliac_vein_R26,
            E = Parameters_Venous1.E_common_iliac_vein_R26,
            r = Parameters_Venous1.r_common_iliac_vein_R26)
        annotation (Placement(transformation(extent={{288,-5},{308,0}})));
          Systemic_vein external_iliac_vein_R28(
            l=Parameters_Venous1.l_external_iliac_vein_R28,
            E=Parameters_Venous1.E_external_iliac_vein_R28,
            r=Parameters_Venous1.r_external_iliac_vein_R28)
            annotation (Placement(transformation(extent={{263,-5},{283,0}})));
        Systemic_vein internal_iliac_vein_T1_R30(
            l = Parameters_Venous1.l_internal_iliac_vein_T1_R30,
            E = Parameters_Venous1.E_internal_iliac_vein_T1_R30,
            r = Parameters_Venous1.r_internal_iliac_vein_T1_R30)
        annotation (Placement(transformation(extent={{60,27},{80,32}})));
        Systemic_vein external_iliac_vein_R32(
            l = Parameters_Venous1.l_external_iliac_vein_R32,
            E = Parameters_Venous1.E_external_iliac_vein_R32,
            r = Parameters_Venous1.r_external_iliac_vein_R32)
        annotation (Placement(transformation(extent={{237,-5},{257,0}})));
        Systemic_vein femoral_vein_R34(
            l = Parameters_Venous1.l_femoral_vein_R34,
            E = Parameters_Venous1.E_femoral_vein_R34,
            r = Parameters_Venous1.r_femoral_vein_R34)
        annotation (Placement(transformation(extent={{212,-5},{232,0}})));
          Systemic_vein femoral_vein_R38(
            l=Parameters_Venous1.l_femoral_vein_R38,
            E=Parameters_Venous1.E_femoral_vein_R38,
            r=Parameters_Venous1.r_femoral_vein_R38)
            annotation (Placement(transformation(extent={{189,-5},{209,0}})));
        Systemic_vein profunda_femoris_vein_T2_R40(
            l = Parameters_Venous1.l_profunda_femoris_vein_T2_R40,
            E = Parameters_Venous1.E_profunda_femoris_vein_T2_R40,
            r = Parameters_Venous1.r_profunda_femoris_vein_T2_R40)
        annotation (Placement(transformation(extent={{60,17},{80,22}})));
        Systemic_vein femoral_vein_R42(
            l = Parameters_Venous1.l_femoral_vein_R42,
            E = Parameters_Venous1.E_femoral_vein_R42,
            r = Parameters_Venous1.r_femoral_vein_R42)
        annotation (Placement(transformation(extent={{163,-5},{183,0}})));
        Systemic_vein femoral_vein_R46(
            l = Parameters_Venous1.l_femoral_vein_R46,
            E = Parameters_Venous1.E_femoral_vein_R46,
            r = Parameters_Venous1.r_femoral_vein_R46)
        annotation (Placement(transformation(extent={{138,-5},{158,0}})));
          Systemic_vein popliteal_vein_R48(
            l=Parameters_Venous1.l_popliteal_vein_R48,
            E=Parameters_Venous1.E_popliteal_vein_R48,
            r=Parameters_Venous1.r_popliteal_vein_R48)
            annotation (Placement(transformation(extent={{113,-5},{133,0}})));
        Systemic_vein anterior_tibial_vein_T4_R50(
            l = Parameters_Venous1.l_anterior_tibial_vein_T4_R50,
            E = Parameters_Venous1.E_anterior_tibial_vein_T4_R50,
            r = Parameters_Venous1.r_anterior_tibial_vein_T4_R50)
        annotation (Placement(transformation(extent={{60,7},{80,12}})));
        Systemic_vein popliteal_vein_R52(
            l = Parameters_Venous1.l_popliteal_vein_R52,
            E = Parameters_Venous1.E_popliteal_vein_R52,
            r = Parameters_Venous1.r_popliteal_vein_R52)
        annotation (Placement(transformation(extent={{87,-5},{107,0}})));
        Systemic_vein posterior_tibial_vein_T6_R54(
            l = Parameters_Venous1.l_posterior_tibial_vein_T6_R54,
            E = Parameters_Venous1.E_posterior_tibial_vein_T6_R54,
            r = Parameters_Venous1.r_posterior_tibial_vein_T6_R54)
        annotation (Placement(transformation(extent={{60,-5},{80,0}})));
          Systemic_vein external_iliac_vein_L58(
            l=Parameters_Venous1.l_external_iliac_vein_L58,
            E=Parameters_Venous1.E_external_iliac_vein_L58,
            r=Parameters_Venous1.r_external_iliac_vein_L58)
            annotation (Placement(transformation(extent={{263,-23},{283,-18}})));
        Systemic_vein internal_iliac_vein_T1_L60(
            l = Parameters_Venous1.l_internal_iliac_vein_T1_L60,
            E = Parameters_Venous1.E_internal_iliac_vein_T1_L60,
            r = Parameters_Venous1.r_internal_iliac_vein_T1_L60)
        annotation (Placement(transformation(extent={{60,-55},{80,-50}})));
        Systemic_vein external_iliac_vein_L62(
            l = Parameters_Venous1.l_external_iliac_vein_L62,
            E = Parameters_Venous1.E_external_iliac_vein_L62,
            r = Parameters_Venous1.r_external_iliac_vein_L62)
        annotation (Placement(transformation(extent={{237,-23},{257,-18}})));
        Systemic_vein femoral_vein_L64(
            l = Parameters_Venous1.l_femoral_vein_L64,
            E = Parameters_Venous1.E_femoral_vein_L64,
            r = Parameters_Venous1.r_femoral_vein_L64)
        annotation (Placement(transformation(extent={{212,-23},{232,-18}})));
          Systemic_vein femoral_vein_L68(
            l=Parameters_Venous1.l_femoral_vein_L68,
            E=Parameters_Venous1.E_femoral_vein_L68,
            r=Parameters_Venous1.r_femoral_vein_L68)
            annotation (Placement(transformation(extent={{189,-23},{209,-18}})));
        Systemic_vein profunda_femoris_vein_T2_L70(
            l = Parameters_Venous1.l_profunda_femoris_vein_T2_L70,
            E = Parameters_Venous1.E_profunda_femoris_vein_T2_L70,
            r = Parameters_Venous1.r_profunda_femoris_vein_T2_L70)
        annotation (Placement(transformation(extent={{60,-43},{80,-38}})));
        Systemic_vein femoral_vein_L72(
            l = Parameters_Venous1.l_femoral_vein_L72,
            E = Parameters_Venous1.E_femoral_vein_L72,
            r = Parameters_Venous1.r_femoral_vein_L72)
        annotation (Placement(transformation(extent={{163,-23},{183,-18}})));
        Systemic_vein femoral_vein_L76(
            l = Parameters_Venous1.l_femoral_vein_L76,
            E = Parameters_Venous1.E_femoral_vein_L76,
            r = Parameters_Venous1.r_femoral_vein_L76)
        annotation (Placement(transformation(extent={{138,-23},{158,-18}})));
          Systemic_vein popliteal_vein_L78(
            l=Parameters_Venous1.l_popliteal_vein_L78,
            E=Parameters_Venous1.E_popliteal_vein_L78,
            r=Parameters_Venous1.r_popliteal_vein_L78)
            annotation (Placement(transformation(extent={{113,-23},{133,-18}})));
        Systemic_vein anterior_tibial_vein_T4_L80(
            l = Parameters_Venous1.l_anterior_tibial_vein_T4_L80,
            E = Parameters_Venous1.E_anterior_tibial_vein_T4_L80,
            r = Parameters_Venous1.r_anterior_tibial_vein_T4_L80)
        annotation (Placement(transformation(extent={{60,-33},{80,-28}})));
        Systemic_vein popliteal_vein_L82(
            l = Parameters_Venous1.l_popliteal_vein_L82,
            E = Parameters_Venous1.E_popliteal_vein_L82,
            r = Parameters_Venous1.r_popliteal_vein_L82)
        annotation (Placement(transformation(extent={{87,-23},{107,-18}})));
        Systemic_vein posterior_tibial_vein_T6_L84(
            l = Parameters_Venous1.l_posterior_tibial_vein_T6_L84,
            E = Parameters_Venous1.E_posterior_tibial_vein_T6_L84,
            r = Parameters_Venous1.r_posterior_tibial_vein_T6_L84)
        annotation (Placement(transformation(extent={{60,-23},{80,-18}})));
          Systemic_vein brachiocephalic_vein_R90(
            l=Parameters_Venous1.l_brachiocephalic_vein_R90,
            E=Parameters_Venous1.E_brachiocephalic_vein_R90,
            r=Parameters_Venous1.r_brachiocephalic_vein_R90)
            annotation (Placement(transformation(extent={{247,163},{267,168}})));
          Systemic_vein brachiocephalic_vein_L124(
            l=Parameters_Venous1.l_brachiocephalic_vein_L124,
            E=Parameters_Venous1.E_brachiocephalic_vein_L124,
            r=Parameters_Venous1.r_brachiocephalic_vein_L124)
            annotation (Placement(transformation(extent={{248,109},{268,114}})));
        Systemic_vein vertebral_vein_R92(
            l = Parameters_Venous1.l_vertebral_vein_R92,
            E = Parameters_Venous1.E_vertebral_vein_R92,
            r = Parameters_Venous1.r_vertebral_vein_R92)
        annotation (Placement(transformation(extent={{61,143},{81,148}})));
          Systemic_vein brachiocephalic_vein_R94(
            l=Parameters_Venous1.l_brachiocephalic_vein_R94,
            E=Parameters_Venous1.E_brachiocephalic_vein_R94,
            r=Parameters_Venous1.r_brachiocephalic_vein_R94)
            annotation (Placement(transformation(extent={{222,163},{242,168}})));
          Systemic_vein subclavian_vein_R96(
            l=Parameters_Venous1.l_subclavian_vein_R96,
            E=Parameters_Venous1.E_subclavian_vein_R96,
            r=Parameters_Venous1.r_subclavian_vein_R96)
            annotation (Placement(transformation(extent={{197,163},{217,168}})));
        Systemic_vein internal_jugular_vein_R122(
            l = Parameters_Venous1.l_internal_jugular_vein_R122,
            E = Parameters_Venous1.E_internal_jugular_vein_R122,
            r = Parameters_Venous1.r_internal_jugular_vein_R122)
        annotation (Placement(transformation(extent={{61,187},{81,192}})));
        Systemic_vein external_jugular_vein_R98(
            l = Parameters_Venous1.l_external_jugular_vein_R98,
            E = Parameters_Venous1.E_external_jugular_vein_R98,
            r = Parameters_Venous1.r_external_jugular_vein_R98)
        annotation (Placement(transformation(extent={{61,175},{81,180}})));
        Systemic_vein subclavian_vein_R100(
            l = Parameters_Venous1.l_subclavian_vein_R100,
            E = Parameters_Venous1.E_subclavian_vein_R100,
            r = Parameters_Venous1.r_subclavian_vein_R100)
        annotation (Placement(transformation(extent={{168,163},{188,168}})));
          Systemic_vein axillary_vein_R102(
            l=Parameters_Venous1.l_axillary_vein_R102,
            E=Parameters_Venous1.E_axillary_vein_R102,
            r=Parameters_Venous1.r_axillary_vein_R102)
            annotation (Placement(transformation(extent={{141,163},{161,168}})));
        Systemic_vein brachial_vein_R104(
            l = Parameters_Venous1.l_brachial_vein_R104,
            E = Parameters_Venous1.E_brachial_vein_R104,
            r = Parameters_Venous1.r_brachial_vein_R104)
        annotation (Placement(transformation(extent={{114,163},{134,168}})));
        Systemic_vein brachial_vein_R114(
            l = Parameters_Venous1.l_brachial_vein_R114,
            E = Parameters_Venous1.E_brachial_vein_R114,
            r = Parameters_Venous1.r_brachial_vein_R114)
        annotation (Placement(transformation(extent={{113,153},{133,158}})));
        Systemic_vein brachial_vein_R108(
            l = Parameters_Venous1.l_brachial_vein_R108,
            E = Parameters_Venous1.E_brachial_vein_R108,
            r = Parameters_Venous1.r_brachial_vein_R108)
        annotation (Placement(transformation(extent={{88,163},{108,168}})));
        Systemic_vein ulnar_vein_T7_R110(
            l = Parameters_Venous1.l_ulnar_vein_T7_R110,
            E = Parameters_Venous1.E_ulnar_vein_T7_R110,
            r = Parameters_Venous1.r_ulnar_vein_T7_R110)
        annotation (Placement(transformation(extent={{61,163},{81,168}})));
        Systemic_vein brachial_vein_R118(
            l = Parameters_Venous1.l_brachial_vein_R118,
            E = Parameters_Venous1.E_brachial_vein_R118,
            r = Parameters_Venous1.r_brachial_vein_R118)
        annotation (Placement(transformation(extent={{88,153},{108,158}})));
        Systemic_vein radial_vein_T3_R120(
            l = Parameters_Venous1.l_radial_vein_T3_R120,
            E = Parameters_Venous1.E_radial_vein_T3_R120,
            r = Parameters_Venous1.r_radial_vein_T3_R120)
        annotation (Placement(transformation(extent={{61,153},{81,158}})));
        Systemic_vein vertebral_vein_L126(
            l = Parameters_Venous1.l_vertebral_vein_L126,
            E = Parameters_Venous1.E_vertebral_vein_L126,
            r = Parameters_Venous1.r_vertebral_vein_L126)
        annotation (Placement(transformation(extent={{60,87},{80,92}})));
          Systemic_vein brachiocephalic_vein_L128(
            l=Parameters_Venous1.l_brachiocephalic_vein_L128,
            E=Parameters_Venous1.E_brachiocephalic_vein_L128,
            r=Parameters_Venous1.r_brachiocephalic_vein_L128)
            annotation (Placement(transformation(extent={{219,109},{239,114}})));
          Systemic_vein subclavian_vein_L130(
            l=Parameters_Venous1.l_subclavian_vein_L130,
            E=Parameters_Venous1.E_subclavian_vein_L130,
            r=Parameters_Venous1.r_subclavian_vein_L130)
            annotation (Placement(transformation(extent={{192,109},{212,114}})));
        Systemic_vein internal_jugular_vein_L156(
            l = Parameters_Venous1.l_internal_jugular_vein_L156,
            E = Parameters_Venous1.E_internal_jugular_vein_L156,
            r = Parameters_Venous1.r_internal_jugular_vein_L156)
        annotation (Placement(transformation(extent={{61,127},{81,132}})));
        Systemic_vein external_jugular_vein_L132(
            l = Parameters_Venous1.l_external_jugular_vein_L132,
            E = Parameters_Venous1.E_external_jugular_vein_L132,
            r = Parameters_Venous1.r_external_jugular_vein_L132)
        annotation (Placement(transformation(extent={{61,119},{81,124}})));
        Systemic_vein subclavian_vein_L134(
            l = Parameters_Venous1.l_subclavian_vein_L134,
            E = Parameters_Venous1.E_subclavian_vein_L134,
            r = Parameters_Venous1.r_subclavian_vein_L134)
        annotation (Placement(transformation(extent={{167,109},{187,114}})));
          Systemic_vein axillary_vein_L136(
            l=Parameters_Venous1.l_axillary_vein_L136,
            E=Parameters_Venous1.E_axillary_vein_L136,
            r=Parameters_Venous1.r_axillary_vein_L136)
            annotation (Placement(transformation(extent={{142,109},{162,114}})));
        Systemic_vein brachial_vein_L138(
            l = Parameters_Venous1.l_brachial_vein_L138,
            E = Parameters_Venous1.E_brachial_vein_L138,
            r = Parameters_Venous1.r_brachial_vein_L138)
        annotation (Placement(transformation(extent={{115,109},{135,114}})));
        Systemic_vein brachial_vein_L148(
            l = Parameters_Venous1.l_brachial_vein_L148,
            E = Parameters_Venous1.E_brachial_vein_L148,
            r = Parameters_Venous1.r_brachial_vein_L148)
        annotation (Placement(transformation(extent={{114,97},{134,102}})));
        Systemic_vein brachial_vein_L142(
            l = Parameters_Venous1.l_brachial_vein_L142,
            E = Parameters_Venous1.E_brachial_vein_L142,
            r = Parameters_Venous1.r_brachial_vein_L142)
        annotation (Placement(transformation(extent={{88,109},{108,114}})));
        Systemic_vein ulnar_vein_T7_L144(
            l = Parameters_Venous1.l_ulnar_vein_T7_L144,
            E = Parameters_Venous1.E_ulnar_vein_T7_L144,
            r = Parameters_Venous1.r_ulnar_vein_T7_L144)
        annotation (Placement(transformation(extent={{61,109},{81,114}})));
        Systemic_vein brachial_vein_L152(
            l = Parameters_Venous1.l_brachial_vein_L152,
            E = Parameters_Venous1.E_brachial_vein_L152,
            r = Parameters_Venous1.r_brachial_vein_L152)
        annotation (Placement(transformation(extent={{88,97},{108,102}})));
        Systemic_vein radial_vein_T3_L154(
            l = Parameters_Venous1.l_radial_vein_T3_L154,
            E = Parameters_Venous1.E_radial_vein_T3_L154,
            r = Parameters_Venous1.r_radial_vein_T3_L154)
        annotation (Placement(transformation(extent={{61,97},{81,102}})));
        ADAN_main.Components.AdanVenousRed._b580e.Parameters_Venous_cellml.Parameters_Venous
          Parameters_Venous1
          annotation (Placement(transformation(extent={{-69,-87},{-49,-82}})));
        Real u_root(unit = "Pa") = ascending_aorta_A.u_in;
        Real v_sup_venacava(unit = "m3.s-1")= superior_vena_cava_C2.v_out;
        Real v_inf_venacava(unit = "m3.s-1")= inferior_vena_cava_C8.v_out;

        replaceable SystemicTissueParameters.SystemicTissueParameters_Calculated
          tissueParameters(
          tissue_pressure=settings.tissues_nominal_pressure,
          arterioles_pressure=settings.tissues_nominal_arterioles_pressure,
          venules_pressure=settings.tissues_nominal_venules_pressure,
          total_zpv=settings.tissues_nominal_zpv,
          stressed_volume=settings.tissues_nominal_stressed_volume,
          cardiac_output=settings.tissues_nominal_cardiac_output)
                           constrainedby
          SystemicTissueParameters.SystemicTissueParameters annotation (
            Placement(transformation(extent={{-140,-98},{-120,-78}})),
            __Dymola_choicesAllMatching=true);
        Systemic_artery                  mesenteric_artery(
          l(displayUnit="cm") = 0.108,
          E(displayUnit="Pa") = 4.00E+05,
          r=3.73E-03) annotation (Placement(transformation(extent={{9,47},{29,52}})));
        Systemic_tissue                          splachnic_tissue(
          Ra=tissueParameters.Ra_splachnic_tissue,
          Rv=tissueParameters.Rv_splachnic_tissue,
          I=tissueParameters.I_splachnic_tissue,
          C=tissueParameters.C_splachnic_tissue,
          zpv=tissueParameters.Zpv_splachnic_tissue,
          nominal_pressure=tissueParameters.tissue_pressure)
          annotation (Placement(transformation(extent={{35,47},{55,52}})));
        Systemic_vein                    splachnic_vein(
          l=1.00E-01,
          r=7.50E-03) annotation (Placement(transformation(extent={{60,47},{80,52}})));
        Systemic_artery_thoracic                  coronary_arteries(
          C(displayUnit="m3/Pa") = 3e-11,
          R(displayUnit="(Pa.s)/m3") = 7e8)
          annotation (Placement(transformation(extent={{10,74},{30,78}})));
        Systemic_tissue                          cardiac_tissue(
          UseOuter_thoracic_pressure=true,
          I=tissueParameters.I_cardiac_tissue,
          C=tissueParameters.C_cardiac_tissue,
          Ra=tissueParameters.Ra_cardiac_tissue,
          Rv=tissueParameters.Rv_cardiac_tissue,
          zpv=tissueParameters.Zpv_cardiac_tissue,
          nominal_pressure=tissueParameters.tissue_pressure)
          annotation (Placement(transformation(extent={{36,74},{56,78}})));
        replaceable
        Systemic_vein coronary_veins(
          UseOuter_thoracic_pressure=true,
          C(displayUnit="m3/Pa") = 7e-10,
          R(displayUnit="(Pa.s)/m3") = 2e8)
          annotation (Placement(transformation(extent={{60,74},{80,78}})));
      equation

        connect(internal_iliac_T1_R218.port_b,internal_iliac_vein_T1_R30.port_a) annotation (Line(points={{55,29.5},
                {60,29.5}},                                                                                                              thickness=1,
            color={28,108,200}));
        connect(internal_iliac_T1_L196.port_b,internal_iliac_vein_T1_L60.port_a) annotation (Line(points={{55,
                -52.5},{60,-52.5}},                                                                                                      thickness=1,
            color={28,108,200}));
        connect(profundus_T2_R224.port_b,profunda_femoris_vein_T2_R40.port_a) annotation (Line(points={{55,19.5},
                {60,19.5}},                                                                                                           thickness=1,
            color={28,108,200}));
        connect(profundus_T2_L202.port_b,profunda_femoris_vein_T2_L70.port_a) annotation (Line(points={{55,
                -40.5},{60,-40.5}},                                                                                                   thickness=1,
            color={28,108,200}));
        connect(anterior_tibial_T3_L208.port_b,anterior_tibial_vein_T4_L80.port_a) annotation (Line(points={{55,
                -30.5},{60,-30.5}},                                                                                                        thickness=1,
            color={28,108,200}));
        connect(anterior_tibial_T3_R230.port_b,anterior_tibial_vein_T4_R50.port_a) annotation (Line(points={{55,9.5},
                {60,9.5}},                                                                                                                 thickness=1,
            color={28,108,200}));
        connect(posterior_tibial_T4_L214.port_b,posterior_tibial_vein_T6_L84.port_a) annotation (Line(points={{55,
                -20.5},{60,-20.5}},                                                                                                          thickness=1,
            color={28,108,200}));
        connect(posterior_tibial_T4_R236.port_b,posterior_tibial_vein_T6_R54.port_a) annotation (Line(points={{55,-2.5},
                {60,-2.5}},                                                                                                                  thickness=1,
            color={28,108,200}));
        connect(radial_T1_R44.port_b,radial_vein_T3_R120.port_a) annotation (Line(points={{55,
                155.5},{61,155.5}},                                                                                      thickness=1,
            color={28,108,200}));
        connect(radial_T1_L92.port_b,radial_vein_T3_L154.port_a) annotation (Line(points={{55,99.5},
                {61,99.5}},                                                                                              thickness=1,
            color={28,108,200}));
        connect(ulnar_T2_R42.port_b,ulnar_vein_T7_R110.port_a) annotation (Line(points={{55,
                165.5},{55,166},{60,166},{60,165.5},{61,165.5}},                                                       thickness=1,
            color={28,108,200}));
        connect(ulnar_T2_L90.port_b,ulnar_vein_T7_L144.port_a) annotation (Line(points={{55,
                111.5},{55,112},{56,112},{56,111.5},{61,111.5}},                                                       thickness=1,
            color={28,108,200}));
        connect(vertebral_R272.port_b,vertebral_vein_R92.port_a) annotation (Line(points={{55,
                145.5},{61,145.5}},                                                                                      thickness=1,
            color={28,108,200}));
        connect(vertebral_L2.port_b,vertebral_vein_L126.port_a) annotation (Line(points={{55,89.5},
                {60,89.5}},                                                                                             thickness=1,
            color={28,108,200}));
        connect(internal_carotid_R8_C.port_b,internal_jugular_vein_R122.port_a) annotation (Line(points={{55,
                189.5},{61,189.5}},                                                                                                     thickness=1,
            color={28,108,200}));
        connect(external_carotid_T2_R26.port_b,external_jugular_vein_R98.port_a) annotation (Line(points={{55,
                177.5},{61,177.5}},                                                                                                      thickness=1,
            color={28,108,200}));
        connect(internal_carotid_L50_C.port_b,internal_jugular_vein_L156.port_a) annotation (Line(points={{55,
                129.5},{61,129.5}},                                                                                                      thickness=1,
            color={28,108,200}));
        connect(external_carotid_T2_L62.port_b,external_jugular_vein_L132.port_a) annotation (Line(points={{55,
                121.5},{61,121.5}},                                                                                                       thickness=1,
            color={28,108,200}));
        connect(celiac_trunk_C116.port_b,hepatic_vein_T1_C10.port_a) annotation (Line(points={{55,59.5},
                {56,59.5},{56,60},{60,60},{60,59.5}},                                                                        thickness=1,
            color={28,108,200}));
        connect(renal_L166.port_b,renal_vein_T1_L22.port_a) annotation (Line(points={{55,
                -64.5},{60,-64.5}},                                                                                 thickness=1,
            color={28,108,200}));
        connect(renal_R178.port_b,renal_vein_T1_R18.port_a) annotation (Line(points={{55,39.5},
                {60,39.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(brachiocephalic_trunk_C4.port_a,aortic_arch_C2.port_b) annotation (Line(points={{-132,
                165.5},{-132,166},{-140,166},{-140,87.5},{-143,87.5}},                                                         thickness=1,
            color={238,46,47}));
        connect(aortic_arch_C46.port_a,aortic_arch_C2.port_b) annotation (Line(points={{-137,
                87.5},{-143,87.5}},                                                                                   thickness=1));
        connect(common_carotid_R6_A.port_a,brachiocephalic_trunk_C4.port_b) annotation (Line(points={{-88,
                189.5},{-112,189.5},{-112,165.5}},                                                                                  thickness=1,
            color={238,46,47}));
        connect(subclavian_R28.port_a,brachiocephalic_trunk_C4.port_b) annotation (Line(points={{-88,
                165.5},{-112,165.5}},                                                                                          thickness=1,
            color={238,46,47}));
        connect(aortic_arch_C64.port_a,aortic_arch_C46.port_b) annotation (Line(points={{-112,
                87.5},{-117,87.5}},                                                                                    thickness=1));
        connect(common_carotid_L48_A.port_a,aortic_arch_C46.port_b) annotation (Line(points={{-111,
                121.5},{-117,121.5},{-117,87.5}},                                                                           thickness=1,
            color={238,46,47}));
        connect(aortic_arch_C94.port_a,aortic_arch_C64.port_b) annotation (Line(points={{-91,
                67.5},{-82,67.5},{-82,88},{-92,88},{-92,87.5}},                                                        thickness=1));
        connect(subclavian_L66.port_a,aortic_arch_C64.port_b) annotation (Line(points={{-87,
                111.5},{-90,111.5},{-90,87.5},{-92,87.5}},                                                            thickness=1,
            color={238,46,47}));
        connect(thoracic_aorta_C108.port_a,thoracic_aorta_C104.port_b) annotation (Line(points={{-193,
                67.5},{-188,67.5}},                                                                                            thickness=1));
        connect(thoracic_aorta_C112.port_a,thoracic_aorta_C108.port_b) annotation (Line(points={{-218,
                67.5},{-213,67.5}},                                                                                            thickness=1));
        connect(abdominal_aorta_C136.port_a,abdominal_aorta_C114.port_b) annotation (Line(points={{-280,
                -2.5},{-287,-2.5}},                                                                                              thickness=1,
            color={238,46,47}));
        connect(celiac_trunk_C116.port_a,abdominal_aorta_C114.port_b) annotation (Line(points={{35,59.5},
                {-284,59.5},{-284,-2.5},{-287,-2.5}},                                                                         thickness=1,
            color={238,46,47}));
        connect(abdominal_aorta_C164.port_a,abdominal_aorta_C136.port_b) annotation (Line(points={{-255,
                -2.5},{-260,-2.5}},                                                                                              thickness=1,
            color={238,46,47}));
        connect(abdominal_aorta_C176.port_a,abdominal_aorta_C164.port_b) annotation (Line(points={{-230,
                -2.5},{-234,-2.5},{-234,-2},{-235,-2.5}},                                                                        thickness=1,
            color={238,46,47}));
        connect(renal_L166.port_a,abdominal_aorta_C164.port_b) annotation (Line(points={{35,
                -64.5},{35,-64},{-232,-64},{-232,-2.5},{-235,-2.5}},                                                   thickness=1,
            color={238,46,47}));
        connect(abdominal_aorta_C188.port_a,abdominal_aorta_C176.port_b) annotation (Line(points={{-203,
                -2.5},{-210,-2.5}},                                                                                              thickness=1,
            color={238,46,47}));
        connect(renal_R178.port_a,abdominal_aorta_C176.port_b) annotation (Line(points={{35,39.5},
                {-210,39.5},{-210,-2.5}},                                                                              thickness=1,
            color={238,46,47}));
        connect(abdominal_aorta_C192.port_a,abdominal_aorta_C188.port_b) annotation (Line(points={{-178,
                -2.5},{-183,-2.5}},                                                                                              thickness=1,
            color={238,46,47}));
        connect(common_iliac_R216.port_a,abdominal_aorta_C192.port_b) annotation (Line(points={{-148,
                -2.5},{-158,-2.5}},                                                                                           thickness=1,
            color={238,46,47}));
        connect(common_iliac_L194.port_a,abdominal_aorta_C192.port_b) annotation (Line(points={{-147,
                -20.5},{-152,-20.5},{-152,-2},{-158,-2},{-158,-2.5}},                                                         thickness=1,
            color={238,46,47}));
        connect(internal_iliac_T1_R218.port_a,common_iliac_R216.port_b) annotation (Line(points={{35,29.5},
                {-128,29.5},{-128,-2.5}},                                                                                       thickness=1,
            color={238,46,47}));
        connect(external_iliac_R220.port_a,common_iliac_R216.port_b) annotation (Line(points={{-120,
                -2.5},{-128,-2.5}},                                                                                          thickness=1,
            color={238,46,47}));
        connect(profundus_T2_R224.port_a,femoral_R222.port_b) annotation (Line(points={{35,19.5},
                {35,20},{-74,20},{-74,-2.5},{-75,-2.5}},                                                              thickness=1,
            color={238,46,47}));
        connect(femoral_R226.port_a,femoral_R222.port_b) annotation (Line(points={{-67,
                -2.5},{-75,-2.5}},                                                                               thickness=1,
            color={238,46,47}));
        connect(anterior_tibial_T3_R230.port_a,popliteal_R228.port_b) annotation (Line(points={{35,9.5},
                {-22,9.5},{-22,-2.5}},                                                                                        thickness=1,
            color={238,46,47}));
        connect(popliteal_R232.port_a,popliteal_R228.port_b) annotation (Line(points={{-16,
                -2.5},{-22,-2.5}},                                                                                   thickness=1,
            color={238,46,47}));
        connect(internal_iliac_T1_L196.port_a,common_iliac_L194.port_b) annotation (Line(points={{35,
                -52.5},{35,-52},{-126,-52},{-126,-20.5},{-127,-20.5}},                                                          thickness=1,
            color={238,46,47}));
        connect(external_iliac_L198.port_a,common_iliac_L194.port_b) annotation (Line(points={{-121,
                -20.5},{-127,-20.5}},                                                                                        thickness=1,
            color={238,46,47}));
        connect(profundus_T2_L202.port_a,femoral_L200.port_b) annotation (Line(points={{35,
                -40.5},{-74,-40.5},{-74,-20.5}},                                                                      thickness=1,
            color={238,46,47}));
        connect(femoral_L204.port_a,femoral_L200.port_b) annotation (Line(points={{-68,
                -20.5},{-74,-20.5}},                                                                             thickness=1,
            color={238,46,47}));
        connect(anterior_tibial_T3_L208.port_a,popliteal_L206.port_b) annotation (Line(points={{35,
                -30.5},{35,-30},{-20,-30},{-20,-20.5},{-21,-20.5}},                                                           thickness=1,
            color={238,46,47}));
        connect(popliteal_L210.port_a,popliteal_L206.port_b) annotation (Line(points={{-16,
                -20.5},{-21,-20.5}},                                                                                 thickness=1,
            color={238,46,47}));
        connect(subclavian_R30.port_a,subclavian_R28.port_b) annotation (Line(points={{-63,
                165.5},{-68,165.5}},                                                                                 thickness=1,
            color={238,46,47}));
        connect(vertebral_R272.port_a,subclavian_R28.port_b) annotation (Line(points={{35,
                145.5},{-68,145.5},{-68,165.5}},                                                                     thickness=1,
            color={238,46,47}));
        connect(ulnar_T2_R36.port_a,brachial_R34.port_b) annotation (Line(points={{10,
                165.5},{3,165.5}},                                                                               thickness=1,
            color={238,46,47}));
        connect(radial_T1_R44.port_a,brachial_R34.port_b) annotation (Line(points={{35,
                155.5},{6,155.5},{6,165.5},{3,165.5}},                                                            thickness=1,
            color={238,46,47}));
        connect(ulnar_T2_R42.port_a,ulnar_T2_R36.port_b) annotation (Line(points={{35,
                165.5},{30,165.5}},                                                                              thickness=1,
            color={238,46,47}));
        connect(subclavian_L78.port_a,subclavian_L66.port_b) annotation (Line(points={{-62,
                111.5},{-67,111.5}},                                                                                 thickness=1,
            color={238,46,47}));
        connect(vertebral_L2.port_a,subclavian_L66.port_b) annotation (Line(points={{35,89.5},
                {-64,89.5},{-64,111.5},{-67,111.5}},                                                               thickness=1,
            color={238,46,47}));
        connect(ulnar_T2_L84.port_a,brachial_L82.port_b) annotation (Line(points={{9,111.5},
                {4,111.5}},                                                                                      thickness=1,
            color={238,46,47}));
        connect(radial_T1_L92.port_a,brachial_L82.port_b) annotation (Line(points={{35,99.5},
                {4,99.5},{4,111.5}},                                                                              thickness=1,
            color={238,46,47}));
        connect(ulnar_T2_L90.port_a,ulnar_T2_L84.port_b) annotation (Line(points={{35,
                111.5},{35,112},{28,112},{28,111.5},{29,111.5}},                                                 thickness=1,
            color={238,46,47}));
        connect(internal_carotid_R8_A.port_a,common_carotid_R6_C.port_b) annotation (Line(points={{-17,
                189.5},{-20,189.5}},                                                                                             thickness=1,
            color={238,46,47}));
        connect(external_carotid_T2_R26.port_a,common_carotid_R6_C.port_b) annotation (Line(points={{35,
                177.5},{-20,177.5},{-20,189.5}},                                                                                   thickness=1,
            color={238,46,47}));
        connect(internal_carotid_L50_A.port_a,common_carotid_L48_D.port_b) annotation (Line(points={{-15,
                129.5},{-18,129.5},{-18,121.5}},                                                                                   thickness=1,
            color={238,46,47}));
        connect(external_carotid_T2_L62.port_a,common_carotid_L48_D.port_b) annotation (Line(points={{35,
                121.5},{-18,121.5}},                                                                                                thickness=1,
            color={238,46,47}));
        connect(ascending_aorta_B.port_a,ascending_aorta_A.port_b) annotation (Line(points={{-238,
                87.5},{-243,87.5}},                                                                                        thickness=1));
        connect(ascending_aorta_C.port_a,ascending_aorta_B.port_b) annotation (Line(points={{-213,
                87.5},{-218,87.5}},                                                                                        thickness=1));
        connect(ascending_aorta_D.port_a,ascending_aorta_C.port_b) annotation (Line(points={{-188,
                87.5},{-193,87.5}},                                                                                        thickness=1));
        connect(aortic_arch_C2.port_a,ascending_aorta_D.port_b) annotation (Line(points={{-163,
                87.5},{-168,87.5}},                                                                                     thickness=1));
        connect(thoracic_aorta_C96.port_a,aortic_arch_C94.port_b) annotation (Line(points={{-118,
                67.5},{-111,67.5}},                                                                                       thickness=1,
            color={0,0,0}));
        connect(abdominal_aorta_C114.port_a,thoracic_aorta_C112.port_b) annotation (Line(points={{-307,
                -2.5},{-307,-2},{-316,-2},{-316,67.5},{-238,67.5}},                                                             thickness=1,
            color={255,0,0}));
        connect(femoral_R222.port_a,external_iliac_R220.port_b) annotation (Line(points={{-95,
                -2.5},{-100,-2.5}},                                                                                     thickness=1,
            color={238,46,47}));
        connect(popliteal_R228.port_a,femoral_R226.port_b) annotation (Line(points={{-42,
                -2.5},{-47,-2.5}},                                                                                 thickness=1,
            color={238,46,47}));
        connect(tibiofibular_trunk_R234.port_a,popliteal_R232.port_b) annotation (Line(points={{9,-2.5},
                {4,-2.5}},                                                                                                    thickness=1,
            color={238,46,47}));
        connect(posterior_tibial_T4_R236.port_a,tibiofibular_trunk_R234.port_b) annotation (Line(points={{35,-2.5},
                {29,-2.5}},                                                                                                             thickness=1,
            color={238,46,47}));
        connect(femoral_L200.port_a,external_iliac_L198.port_b) annotation (Line(points={{-94,
                -20.5},{-101,-20.5}},                                                                                   thickness=1,
            color={238,46,47}));
        connect(popliteal_L206.port_a,femoral_L204.port_b) annotation (Line(points={{-41,
                -20.5},{-48,-20.5}},                                                                               thickness=1,
            color={238,46,47}));
        connect(tibiofibular_trunk_L212.port_a,popliteal_L210.port_b) annotation (Line(points={{9,-20.5},
                {4,-20.5}},                                                                                                   thickness=1,
            color={238,46,47}));
        connect(posterior_tibial_T4_L214.port_a,tibiofibular_trunk_L212.port_b) annotation (Line(points={{35,
                -20.5},{29,-20.5}},                                                                                                     thickness=1,
            color={238,46,47}));
        connect(axillary_R32.port_a,subclavian_R30.port_b) annotation (Line(points={{-40,
                165.5},{-43,165.5}},                                                                               thickness=1,
            color={238,46,47}));
        connect(brachial_R34.port_a,axillary_R32.port_b) annotation (Line(points={{-17,
                165.5},{-20,165.5}},                                                                             thickness=1,
            color={238,46,47}));
        connect(axillary_L80.port_a,subclavian_L78.port_b) annotation (Line(points={{-39,
                111.5},{-42,111.5}},                                                                               thickness=1,
            color={238,46,47}));
        connect(brachial_L82.port_a,axillary_L80.port_b) annotation (Line(points={{-16,
                111.5},{-19,111.5}},                                                                             thickness=1,
            color={238,46,47}));
        connect(common_carotid_R6_B.port_a,common_carotid_R6_A.port_b) annotation (Line(points={{-63,
                189.5},{-68,189.5}},                                                                                           thickness=1,
            color={238,46,47}));
        connect(common_carotid_R6_C.port_a,common_carotid_R6_B.port_b) annotation (Line(points={{-40,
                189.5},{-43,189.5}},                                                                                           thickness=1,
            color={238,46,47}));
        connect(internal_carotid_R8_B.port_a,internal_carotid_R8_A.port_b) annotation (Line(points={{10,
                189.5},{3,189.5}},                                                                                                 thickness=1,
            color={238,46,47}));
        connect(internal_carotid_R8_C.port_a,internal_carotid_R8_B.port_b) annotation (Line(points={{35,
                189.5},{30,189.5}},                                                                                                thickness=1,
            color={238,46,47}));
        connect(common_carotid_L48_B.port_a,common_carotid_L48_A.port_b) annotation (Line(points={{-86,
                121.5},{-91,121.5}},                                                                                             thickness=1,
            color={238,46,47}));
        connect(common_carotid_L48_C.port_a,common_carotid_L48_B.port_b) annotation (Line(points={{-63,
                121.5},{-66,121.5}},                                                                                             thickness=1,
            color={238,46,47}));
        connect(common_carotid_L48_D.port_a,common_carotid_L48_C.port_b) annotation (Line(points={{-38,
                121.5},{-43,121.5}},                                                                                             thickness=1,
            color={238,46,47}));
        connect(internal_carotid_L50_B.port_a,internal_carotid_L50_A.port_b) annotation (Line(points={{10,
                129.5},{5,129.5}},                                                                                                   thickness=1,
            color={238,46,47}));
        connect(internal_carotid_L50_C.port_a,internal_carotid_L50_B.port_b) annotation (Line(points={{35,
                129.5},{30,129.5}},                                                                                                  thickness=1,
            color={238,46,47}));
        connect(superior_vena_cava_C88.port_b,superior_vena_cava_C2.port_a) annotation (Line(points={{292,
                165.5},{297,165.5}},                                                                                                thickness=1,
            color={28,108,200}));
        connect(brachiocephalic_vein_R90.port_b,superior_vena_cava_C88.port_a) annotation (Line(points={{267,
                165.5},{272,165.5}},                                                                                                   thickness=1,
            color={28,108,200}));
        connect(brachiocephalic_vein_L124.port_b,superior_vena_cava_C88.port_a) annotation (Line(points={{268,
                111.5},{268,112},{272,112},{272,165.5}},                                                                                thickness=1,
            color={28,108,200}));
        connect(vertebral_vein_R92.port_b,brachiocephalic_vein_R90.port_a) annotation (Line(points={{81,
                145.5},{247,145.5},{247,165.5}},                                                                                   thickness=1,
            color={28,108,200}));
        connect(brachiocephalic_vein_R94.port_b,brachiocephalic_vein_R90.port_a) annotation (Line(points={{242,
                165.5},{247,165.5}},                                                                                                     thickness=1,
            color={28,108,200}));
        connect(subclavian_vein_R96.port_b,brachiocephalic_vein_R94.port_a) annotation (Line(points={{217,
                165.5},{222,165.5}},                                                                                                thickness=1,
            color={28,108,200}));
        connect(internal_jugular_vein_R122.port_b,brachiocephalic_vein_R94.port_a) annotation (Line(points={{81,
                189.5},{222,189.5},{222,165.5}},                                                                                           thickness=1,
            color={28,108,200}));
        connect(external_jugular_vein_R98.port_b,subclavian_vein_R96.port_a) annotation (Line(points={{81,
                177.5},{196,177.5},{196,165.5},{197,165.5}},                                                                         thickness=1,
            color={28,108,200}));
        connect(subclavian_vein_R100.port_b,subclavian_vein_R96.port_a) annotation (Line(points={{188,
                165.5},{197,165.5}},                                                                                            thickness=1,
            color={28,108,200}));
        connect(brachial_vein_R104.port_b,axillary_vein_R102.port_a) annotation (Line(points={{134,
                165.5},{141,165.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(brachial_vein_R114.port_b,axillary_vein_R102.port_a) annotation (Line(points={{133,
                155.5},{134,155.5},{134,156},{136,156},{136,165.5},{141,165.5}},                                             thickness=1,
            color={28,108,200}));
        connect(ulnar_vein_T7_R110.port_b,brachial_vein_R108.port_a) annotation (Line(points={{81,
                165.5},{88,165.5}},                                                                                          thickness=1,
            color={28,108,200}));
        connect(vertebral_vein_L126.port_b,brachiocephalic_vein_L124.port_a) annotation (Line(points={{80,89.5},
                {80,90},{244,90},{244,111.5},{248,111.5}},                                                                           thickness=1,
            color={28,108,200}));
        connect(brachiocephalic_vein_L128.port_b,brachiocephalic_vein_L124.port_a) annotation (Line(points={{239,
                111.5},{248,111.5}},                                                                                                       thickness=1,
            color={28,108,200}));
        connect(subclavian_vein_L130.port_b,brachiocephalic_vein_L128.port_a) annotation (Line(points={{212,
                111.5},{219,111.5}},                                                                                                  thickness=1,
            color={28,108,200}));
        connect(internal_jugular_vein_L156.port_b,brachiocephalic_vein_L128.port_a) annotation (Line(points={{81,
                129.5},{216,129.5},{216,111.5},{219,111.5}},                                                                                thickness=1,
            color={28,108,200}));
        connect(external_jugular_vein_L132.port_b,subclavian_vein_L130.port_a) annotation (Line(points={{81,
                121.5},{192,121.5},{192,111.5}},                                                                                       thickness=1,
            color={28,108,200}));
        connect(subclavian_vein_L134.port_b,subclavian_vein_L130.port_a) annotation (Line(points={{187,
                111.5},{192,111.5}},                                                                                             thickness=1,
            color={28,108,200}));
        connect(brachial_vein_L138.port_b,axillary_vein_L136.port_a) annotation (Line(points={{135,
                111.5},{142,111.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(brachial_vein_L148.port_b,axillary_vein_L136.port_a) annotation (Line(points={{134,
                99.5},{134,100},{138,100},{138,111.5},{142,111.5}},                                                          thickness=1,
            color={28,108,200}));
        connect(ulnar_vein_T7_L144.port_b,brachial_vein_L142.port_a) annotation (Line(points={{81,
                111.5},{88,111.5}},                                                                                          thickness=1,
            color={28,108,200}));
        connect(hepatic_vein_T1_C10.port_b,inferior_vena_cava_C8.port_a) annotation (Line(points={{80,59.5},
                {406,59.5},{406,60},{408,60},{408,-2.5},{409,-2.5}},                                                             thickness=1,
            color={28,108,200}));
        connect(inferior_vena_cava_C12.port_b,inferior_vena_cava_C8.port_a) annotation (Line(points={{403,
                -2.5},{403,-2},{409,-2},{409,-2.5}},                                                                                thickness=1,
            color={28,108,200}));
        connect(renal_vein_T1_R18.port_b,inferior_vena_cava_C16.port_a) annotation (Line(points={{80,39.5},
                {360,39.5},{360,-2.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(inferior_vena_cava_C20.port_b,inferior_vena_cava_C16.port_a) annotation (Line(points={{356,
                -2.5},{360,-2.5}},                                                                                                   thickness=1,
            color={28,108,200}));
        connect(renal_vein_T1_L22.port_b,inferior_vena_cava_C20.port_a) annotation (Line(points={{80,
                -64.5},{336,-64.5},{336,-2.5}},                                                                                 thickness=1,
            color={28,108,200}));
        connect(inferior_vena_cava_C24.port_b,inferior_vena_cava_C20.port_a) annotation (Line(points={{332,
                -2.5},{336,-2.5}},                                                                                                   thickness=1,
            color={28,108,200}));
        connect(common_iliac_vein_R26.port_b,inferior_vena_cava_C24.port_a) annotation (Line(points={{308,
                -2.5},{312,-2.5}},                                                                                                  thickness=1,
            color={28,108,200}));
        connect(common_iliac_vein_L56.port_b,inferior_vena_cava_C24.port_a) annotation (Line(points={{307,
                -20.5},{312,-20.5},{312,-2.5}},                                                                                     thickness=1,
            color={28,108,200}));
        connect(internal_iliac_vein_T1_R30.port_b,external_iliac_vein_R28.port_a) annotation (Line(points={{80,29.5},
                {80,30},{260,30},{260,-2.5},{263,-2.5}},                                                                                  thickness=1,
            color={28,108,200}));
        connect(external_iliac_vein_R32.port_b,external_iliac_vein_R28.port_a) annotation (Line(points={{257,
                -2.5},{263,-2.5}},                                                                                                     thickness=1,
            color={28,108,200}));
        connect(femoral_vein_R38.port_b,femoral_vein_R34.port_a) annotation (Line(points={{209,
                -2.5},{212,-2.5}},                                                                                       thickness=1,
            color={28,108,200}));
        connect(profunda_femoris_vein_T2_R40.port_b,femoral_vein_R38.port_a) annotation (Line(points={{80,19.5},
                {188,19.5},{188,-2.5},{189,-2.5}},                                                                                   thickness=1,
            color={28,108,200}));
        connect(femoral_vein_R42.port_b,femoral_vein_R38.port_a) annotation (Line(points={{183,
                -2.5},{189,-2.5}},                                                                                       thickness=1,
            color={28,108,200}));
        connect(anterior_tibial_vein_T4_R50.port_b,popliteal_vein_R48.port_a) annotation (Line(points={{80,9.5},
                {112,9.5},{112,-2.5},{113,-2.5}},                                                                                     thickness=1,
            color={28,108,200}));
        connect(popliteal_vein_R52.port_b,popliteal_vein_R48.port_a) annotation (Line(points={{107,
                -2.5},{113,-2.5}},                                                                                           thickness=1,
            color={28,108,200}));
        connect(internal_iliac_vein_T1_L60.port_b,external_iliac_vein_L58.port_a) annotation (Line(points={{80,
                -52.5},{262,-52.5},{262,-20.5},{263,-20.5}},                                                                              thickness=1,
            color={28,108,200}));
        connect(external_iliac_vein_L62.port_b,external_iliac_vein_L58.port_a) annotation (Line(points={{257,
                -20.5},{263,-20.5}},                                                                                                   thickness=1,
            color={28,108,200}));
        connect(femoral_vein_L68.port_b,femoral_vein_L64.port_a) annotation (Line(points={{209,
                -20.5},{212,-20.5}},                                                                                     thickness=1,
            color={28,108,200}));
        connect(profunda_femoris_vein_T2_L70.port_b,femoral_vein_L68.port_a) annotation (Line(points={{80,
                -40.5},{188,-40.5},{188,-20.5},{189,-20.5}},                                                                         thickness=1,
            color={28,108,200}));
        connect(femoral_vein_L72.port_b,femoral_vein_L68.port_a) annotation (Line(points={{183,
                -20.5},{189,-20.5}},                                                                                     thickness=1,
            color={28,108,200}));
        connect(anterior_tibial_vein_T4_L80.port_b,popliteal_vein_L78.port_a) annotation (Line(points={{80,
                -30.5},{112,-30.5},{112,-20.5},{113,-20.5}},                                                                          thickness=1,
            color={28,108,200}));
        connect(popliteal_vein_L82.port_b,popliteal_vein_L78.port_a) annotation (Line(points={{107,
                -20.5},{113,-20.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(brachial_vein_R118.port_b,brachial_vein_R114.port_a) annotation (Line(points={{108,
                155.5},{113,155.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(brachial_vein_R108.port_b,brachial_vein_R104.port_a) annotation (Line(points={{108,
                165.5},{114,165.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(brachial_vein_L142.port_b,brachial_vein_L138.port_a) annotation (Line(points={{108,
                111.5},{115,111.5}},                                                                                         thickness=1,
            color={28,108,200}));
        connect(brachial_vein_L152.port_b,brachial_vein_L148.port_a) annotation (Line(points={{108,
                99.5},{114,99.5}},                                                                                           thickness=1,
            color={28,108,200}));
        connect(inferior_vena_cava_C16.port_b,inferior_vena_cava_C12.port_a) annotation (Line(points={{380,
                -2.5},{383,-2.5}},                                                                                                   thickness=1,
            color={28,108,200}));
        connect(axillary_vein_R102.port_b,subclavian_vein_R100.port_a) annotation (Line(points={{161,
                165.5},{168,165.5}},                                                                                           thickness=1,
            color={28,108,200}));
        connect(radial_vein_T3_R120.port_b,brachial_vein_R118.port_a) annotation (Line(points={{81,
                155.5},{88,155.5}},                                                                                           thickness=1,
            color={28,108,200}));
        connect(axillary_vein_L136.port_b,subclavian_vein_L134.port_a) annotation (Line(points={{162,
                111.5},{167,111.5}},                                                                                           thickness=1,
            color={28,108,200}));
        connect(radial_vein_T3_L154.port_b,brachial_vein_L152.port_a) annotation (Line(points={{81,99.5},
                {88,99.5}},                                                                                                   thickness=1,
            color={28,108,200}));
        connect(external_iliac_vein_R28.port_b,common_iliac_vein_R26.port_a) annotation (Line(points={{283,
                -2.5},{283,-2},{288,-2},{288,-2.5}},                                                                                 thickness=1,
            color={28,108,200}));
        connect(femoral_vein_R34.port_b,external_iliac_vein_R32.port_a) annotation (Line(points={{232,
                -2.5},{237,-2.5}},                                                                                              thickness=1,
            color={28,108,200}));
        connect(popliteal_vein_R48.port_b,femoral_vein_R46.port_a) annotation (Line(points={{133,
                -2.5},{133,-2},{138,-2},{138,-2.5}},                                                                       thickness=1,
            color={28,108,200}));
        connect(posterior_tibial_vein_T6_R54.port_b,popliteal_vein_R52.port_a) annotation (Line(points={{80,-2.5},
                {87,-2.5}},                                                                                                            thickness=1,
            color={28,108,200}));
        connect(external_iliac_vein_L58.port_b,common_iliac_vein_L56.port_a) annotation (Line(points={{283,
                -20.5},{287,-20.5}},                                                                                                 thickness=1,
            color={28,108,200}));
        connect(femoral_vein_L64.port_b,external_iliac_vein_L62.port_a) annotation (Line(points={{232,
                -20.5},{237,-20.5}},                                                                                            thickness=1,
            color={28,108,200}));
        connect(popliteal_vein_L78.port_b,femoral_vein_L76.port_a) annotation (Line(points={{133,
                -20.5},{138,-20.5}},                                                                                       thickness=1,
            color={28,108,200}));
        connect(posterior_tibial_vein_T6_L84.port_b,popliteal_vein_L82.port_a) annotation (Line(points={{80,
                -20.5},{87,-20.5}},                                                                                                    thickness=1,
            color={28,108,200}));
        connect(femoral_vein_R46.port_b,femoral_vein_R42.port_a) annotation (Line(points={{158,
                -2.5},{163,-2.5}},                                                                                       thickness=1,
            color={28,108,200}));
        connect(femoral_vein_L76.port_b,femoral_vein_L72.port_a) annotation (Line(points={{158,
                -20.5},{163,-20.5}},                                                                                     thickness=1,
            color={28,108,200}));
        connect(thoracic_aorta_C96.port_b, thoracic_aorta_C100.port_a) annotation (Line(
            points={{-138,67.5},{-143,67.5}},
            color={0,0,0},
            thickness=1));
        connect(thoracic_aorta_C104.port_a, thoracic_aorta_C100.port_b) annotation (Line(
            points={{-168,67.5},{-163,67.5}},
            color={0,0,0},
            thickness=1));
        connect(mesenteric_artery.port_a, abdominal_aorta_C136.port_b) annotation (
            Line(
            points={{9,49.5},{-258,49.5},{-258,-2.5},{-260,-2.5}},
            color={0,0,0},
            thickness=1));
        connect(splachnic_tissue.port_a, mesenteric_artery.port_b) annotation (Line(
            points={{35,49.5},{29,49.5}},
            color={0,0,0},
            thickness=1));
        connect(splachnic_tissue.port_b, splachnic_vein.port_a) annotation (Line(
            points={{55,49.5},{60,49.5}},
            color={0,0,0},
            thickness=1));
        connect(splachnic_vein.port_b, inferior_vena_cava_C12.port_a) annotation (
            Line(
            points={{80,49.5},{196,49.5},{196,50},{383,50},{383,-2.5}},
            color={0,0,0},
            thickness=1));
        connect(cardiac_tissue.port_b, coronary_veins.port_a) annotation (Line(
            points={{56,76},{60,76}},
            color={0,0,0},
            thickness=1));
        connect(cardiac_tissue.port_a, coronary_arteries.port_b) annotation (Line(
            points={{36,76},{30,76}},
            color={0,0,0},
            thickness=1));
        connect(ascending_aorta_B.port_b, coronary_arteries.port_a) annotation (Line(
            points={{-218,87.5},{-218,76},{10,76}},
            color={0,0,0},
            thickness=1));
        connect(coronary_veins.port_b, inferior_vena_cava_C8.port_b) annotation (Line(
            points={{80,76},{429,76},{429,-2.5}},
            color={0,0,0},
            thickness=1));
          annotation (Diagram(coordinateSystem(extent={{-320,-100},{440,200}})),
              Icon(coordinateSystem(extent={{-320,-100},{440,200}}), graphics={
                        Bitmap(extent={{-320,-100},{440,200}}, fileName=
                    "modelica://Physiolibrary/Resources/Icons/perfusion.png")}));
      end Systemic_base;

      partial model Systemic_volumes
        extends ADAN_main.Components.AdanVenousRed.Systemic_base;
      Physiolibrary.Types.Volume total_volume = volume_arterial + volume_peripheral + volume_venous;
        Physiolibrary.Types.Volume volume_arterial = ascending_aorta_A.volume +
          ascending_aorta_B.volume +
          ascending_aorta_C.volume +
          ascending_aorta_D.volume +
          aortic_arch_C2.volume +
          brachiocephalic_trunk_C4.volume +
          aortic_arch_C46.volume +
          aortic_arch_C64.volume +
          aortic_arch_C94.volume +
          thoracic_aorta_C96.volume +
          thoracic_aorta_C100.volume +
          thoracic_aorta_C104.volume +
          thoracic_aorta_C108.volume +
          thoracic_aorta_C112.volume +
          abdominal_aorta_C114.volume +
          abdominal_aorta_C136.volume +
          abdominal_aorta_C164.volume +
          abdominal_aorta_C176.volume +
          abdominal_aorta_C188.volume +
          abdominal_aorta_C192.volume +
          common_iliac_R216.volume +
          external_iliac_R220.volume +
          femoral_R222.volume +
          femoral_R226.volume +
          popliteal_R228.volume +
          popliteal_R232.volume +
          tibiofibular_trunk_R234.volume +
          common_iliac_L194.volume +
          external_iliac_L198.volume +
          femoral_L200.volume +
          femoral_L204.volume +
          popliteal_L206.volume +
          popliteal_L210.volume +
          tibiofibular_trunk_L212.volume +
          subclavian_R28.volume +
          subclavian_R30.volume +
          axillary_R32.volume +
          brachial_R34.volume +
          ulnar_T2_R36.volume +
          subclavian_L66.volume +
          subclavian_L78.volume +
          axillary_L80.volume +
          brachial_L82.volume +
          ulnar_T2_L84.volume +
          common_carotid_R6_A.volume +
          common_carotid_R6_B.volume +
          common_carotid_R6_C.volume +
          internal_carotid_R8_A.volume +
          internal_carotid_R8_B.volume +
          common_carotid_L48_A.volume +
          common_carotid_L48_B.volume +
          common_carotid_L48_C.volume +
          common_carotid_L48_D.volume +
          internal_carotid_L50_A.volume +
          internal_carotid_L50_B.volume +
          coronary_arteries.volume +
          mesenteric_artery.volume;
        Physiolibrary.Types.Volume volume_peripheral = celiac_trunk_C116.volume +
          renal_L166.volume +
          renal_R178.volume +
          internal_iliac_T1_R218.volume +
          profundus_T2_R224.volume +
          anterior_tibial_T3_R230.volume +
          posterior_tibial_T4_R236.volume +
          internal_iliac_T1_L196.volume +
          profundus_T2_L202.volume +
          anterior_tibial_T3_L208.volume +
          posterior_tibial_T4_L214.volume +
          ulnar_T2_R42.volume +
          radial_T1_R44.volume +
          ulnar_T2_L90.volume +
          radial_T1_L92.volume +
          internal_carotid_R8_C.volume +
          external_carotid_T2_R26.volume +
          internal_carotid_L50_C.volume +
          external_carotid_T2_L62.volume +
          vertebral_L2.volume +
          vertebral_R272.volume +
          cardiac_tissue.volume +
          splachnic_tissue.volume;
        Physiolibrary.Types.Volume volume_venous = superior_vena_cava_C2.volume +
          superior_vena_cava_C88.volume +
          inferior_vena_cava_C8.volume +
          hepatic_vein_T1_C10.volume +
          inferior_vena_cava_C12.volume +
          inferior_vena_cava_C16.volume +
          renal_vein_T1_R18.volume +
          inferior_vena_cava_C20.volume +
          renal_vein_T1_L22.volume +
          inferior_vena_cava_C24.volume +
          common_iliac_vein_L56.volume +
          common_iliac_vein_R26.volume +
          external_iliac_vein_R28.volume +
          internal_iliac_vein_T1_R30.volume +
          external_iliac_vein_R32.volume +
          femoral_vein_R34.volume +
          femoral_vein_R38.volume +
          profunda_femoris_vein_T2_R40.volume +
          femoral_vein_R42.volume +
          femoral_vein_R46.volume +
          popliteal_vein_R48.volume +
          anterior_tibial_vein_T4_R50.volume +
          popliteal_vein_R52.volume +
          posterior_tibial_vein_T6_R54.volume +
          external_iliac_vein_L58.volume +
          internal_iliac_vein_T1_L60.volume +
          external_iliac_vein_L62.volume +
          femoral_vein_L64.volume +
          femoral_vein_L68.volume +
          profunda_femoris_vein_T2_L70.volume +
          femoral_vein_L72.volume +
          femoral_vein_L76.volume +
          popliteal_vein_L78.volume +
          anterior_tibial_vein_T4_L80.volume +
          popliteal_vein_L82.volume +
          posterior_tibial_vein_T6_L84.volume +
          brachiocephalic_vein_R90.volume +
          brachiocephalic_vein_L124.volume +
          vertebral_vein_R92.volume +
          brachiocephalic_vein_R94.volume +
          subclavian_vein_R96.volume +
          internal_jugular_vein_R122.volume +
          external_jugular_vein_R98.volume +
          subclavian_vein_R100.volume +
          axillary_vein_R102.volume +
          brachial_vein_R104.volume +
          brachial_vein_R114.volume +
          brachial_vein_R108.volume +
          ulnar_vein_T7_R110.volume +
          brachial_vein_R118.volume +
          radial_vein_T3_R120.volume +
          vertebral_vein_L126.volume +
          brachiocephalic_vein_L128.volume +
          subclavian_vein_L130.volume +
          internal_jugular_vein_L156.volume +
          external_jugular_vein_L132.volume +
          subclavian_vein_L134.volume +
          axillary_vein_L136.volume +
          brachial_vein_L138.volume +
          brachial_vein_L148.volume +
          brachial_vein_L142.volume +
          ulnar_vein_T7_L144.volume +
          brachial_vein_L152.volume +
          radial_vein_T3_L154.volume +
          coronary_veins.volume +
          splachnic_vein.volume;

        Physiolibrary.Types.Volume zpv_arterial = ascending_aorta_A.zpv +
          ascending_aorta_B.zpv +
          ascending_aorta_C.zpv +
          ascending_aorta_D.zpv +
          aortic_arch_C2.zpv +
          brachiocephalic_trunk_C4.zpv +
          aortic_arch_C46.zpv +
          aortic_arch_C64.zpv +
          aortic_arch_C94.zpv +
          thoracic_aorta_C96.zpv +
          thoracic_aorta_C100.zpv +
          thoracic_aorta_C104.zpv +
          thoracic_aorta_C108.zpv +
          thoracic_aorta_C112.zpv +
          abdominal_aorta_C114.zpv +
          abdominal_aorta_C136.zpv +
          abdominal_aorta_C164.zpv +
          abdominal_aorta_C176.zpv +
          abdominal_aorta_C188.zpv +
          abdominal_aorta_C192.zpv +
          common_iliac_R216.zpv +
          external_iliac_R220.zpv +
          femoral_R222.zpv +
          femoral_R226.zpv +
          popliteal_R228.zpv +
          popliteal_R232.zpv +
          tibiofibular_trunk_R234.zpv +
          common_iliac_L194.zpv +
          external_iliac_L198.zpv +
          femoral_L200.zpv +
          femoral_L204.zpv +
          popliteal_L206.zpv +
          popliteal_L210.zpv +
          tibiofibular_trunk_L212.zpv +
          subclavian_R28.zpv +
          subclavian_R30.zpv +
          axillary_R32.zpv +
          brachial_R34.zpv +
          ulnar_T2_R36.zpv +
          subclavian_L66.zpv +
          subclavian_L78.zpv +
          axillary_L80.zpv +
          brachial_L82.zpv +
          ulnar_T2_L84.zpv +
          common_carotid_R6_A.zpv +
          common_carotid_R6_B.zpv +
          common_carotid_R6_C.zpv +
          internal_carotid_R8_A.zpv +
          internal_carotid_R8_B.zpv +
          common_carotid_L48_A.zpv +
          common_carotid_L48_B.zpv +
          common_carotid_L48_C.zpv +
          common_carotid_L48_D.zpv +
          internal_carotid_L50_A.zpv +
          internal_carotid_L50_B.zpv +
          coronary_arteries.zpv +
          mesenteric_artery.zpv;
        Physiolibrary.Types.Volume zpv_peripheral = celiac_trunk_C116.zpv +
          renal_L166.zpv +
          renal_R178.zpv +
          internal_iliac_T1_R218.zpv +
          profundus_T2_R224.zpv +
          anterior_tibial_T3_R230.zpv +
          posterior_tibial_T4_R236.zpv +
          internal_iliac_T1_L196.zpv +
          profundus_T2_L202.zpv +
          anterior_tibial_T3_L208.zpv +
          posterior_tibial_T4_L214.zpv +
          ulnar_T2_R42.zpv +
          radial_T1_R44.zpv +
          ulnar_T2_L90.zpv +
          radial_T1_L92.zpv +
          internal_carotid_R8_C.zpv +
          external_carotid_T2_R26.zpv +
          internal_carotid_L50_C.zpv +
          external_carotid_T2_L62.zpv +
          vertebral_L2.zpv +
          vertebral_R272.zpv  +
          cardiac_tissue.zpv +
          splachnic_tissue.zpv;
        Physiolibrary.Types.Volume zpv_venous = superior_vena_cava_C2.zpv +
          superior_vena_cava_C88.zpv +
          inferior_vena_cava_C8.zpv +
          hepatic_vein_T1_C10.zpv +
          inferior_vena_cava_C12.zpv +
          inferior_vena_cava_C16.zpv +
          renal_vein_T1_R18.zpv +
          inferior_vena_cava_C20.zpv +
          renal_vein_T1_L22.zpv +
          inferior_vena_cava_C24.zpv +
          common_iliac_vein_L56.zpv +
          common_iliac_vein_R26.zpv +
          external_iliac_vein_R28.zpv +
          internal_iliac_vein_T1_R30.zpv +
          external_iliac_vein_R32.zpv +
          femoral_vein_R34.zpv +
          femoral_vein_R38.zpv +
          profunda_femoris_vein_T2_R40.zpv +
          femoral_vein_R42.zpv +
          femoral_vein_R46.zpv +
          popliteal_vein_R48.zpv +
          anterior_tibial_vein_T4_R50.zpv +
          popliteal_vein_R52.zpv +
          posterior_tibial_vein_T6_R54.zpv +
          external_iliac_vein_L58.zpv +
          internal_iliac_vein_T1_L60.zpv +
          external_iliac_vein_L62.zpv +
          femoral_vein_L64.zpv +
          femoral_vein_L68.zpv +
          profunda_femoris_vein_T2_L70.zpv +
          femoral_vein_L72.zpv +
          femoral_vein_L76.zpv +
          popliteal_vein_L78.zpv +
          anterior_tibial_vein_T4_L80.zpv +
          popliteal_vein_L82.zpv +
          posterior_tibial_vein_T6_L84.zpv +
          brachiocephalic_vein_R90.zpv +
          brachiocephalic_vein_L124.zpv +
          vertebral_vein_R92.zpv +
          brachiocephalic_vein_R94.zpv +
          subclavian_vein_R96.zpv +
          internal_jugular_vein_R122.zpv +
          external_jugular_vein_R98.zpv +
          subclavian_vein_R100.zpv +
          axillary_vein_R102.zpv +
          brachial_vein_R104.zpv +
          brachial_vein_R114.zpv +
          brachial_vein_R108.zpv +
          ulnar_vein_T7_R110.zpv +
          brachial_vein_R118.zpv +
          radial_vein_T3_R120.zpv +
          vertebral_vein_L126.zpv +
          brachiocephalic_vein_L128.zpv +
          subclavian_vein_L130.zpv +
          internal_jugular_vein_L156.zpv +
          external_jugular_vein_L132.zpv +
          subclavian_vein_L134.zpv +
          axillary_vein_L136.zpv +
          brachial_vein_L138.zpv +
          brachial_vein_L148.zpv +
          brachial_vein_L142.zpv +
          ulnar_vein_T7_L144.zpv +
          brachial_vein_L152.zpv +
          radial_vein_T3_L154.zpv +
          coronary_veins.zpv +
          splachnic_vein.zpv;

       // Physiolibrary.Types.Volume zpvs_venous =   zpv_venous * ZPV_effect;

      Modelica.SIunits.Length l_arterial = ascending_aorta_A.l +
          ascending_aorta_B.l +
          ascending_aorta_C.l +
          ascending_aorta_D.l +
          aortic_arch_C2.l +
          brachiocephalic_trunk_C4.l +
          aortic_arch_C46.l +
          aortic_arch_C64.l +
          aortic_arch_C94.l +
          thoracic_aorta_C96.l +
          thoracic_aorta_C100.l +
          thoracic_aorta_C104.l +
          thoracic_aorta_C108.l +
          thoracic_aorta_C112.l +
          abdominal_aorta_C114.l +
          abdominal_aorta_C136.l +
          abdominal_aorta_C164.l +
          abdominal_aorta_C176.l +
          abdominal_aorta_C188.l +
          abdominal_aorta_C192.l +
          common_iliac_R216.l +
          external_iliac_R220.l +
          femoral_R222.l +
          femoral_R226.l +
          popliteal_R228.l +
          popliteal_R232.l +
          tibiofibular_trunk_R234.l +
          common_iliac_L194.l +
          external_iliac_L198.l +
          femoral_L200.l +
          femoral_L204.l +
          popliteal_L206.l +
          popliteal_L210.l +
          tibiofibular_trunk_L212.l +
          subclavian_R28.l +
          subclavian_R30.l +
          axillary_R32.l +
          brachial_R34.l +
          ulnar_T2_R36.l +
          subclavian_L66.l +
          subclavian_L78.l +
          axillary_L80.l +
          brachial_L82.l +
          ulnar_T2_L84.l +
          common_carotid_R6_A.l +
          common_carotid_R6_B.l +
          common_carotid_R6_C.l +
          internal_carotid_R8_A.l +
          internal_carotid_R8_B.l +
          common_carotid_L48_A.l +
          common_carotid_L48_B.l +
          common_carotid_L48_C.l +
          common_carotid_L48_D.l +
          internal_carotid_L50_A.l +
          internal_carotid_L50_B.l +
          coronary_arteries.l +
          mesenteric_artery.l;

      Modelica.SIunits.Length l_venous = superior_vena_cava_C2.l +
          superior_vena_cava_C88.l +
          inferior_vena_cava_C8.l +
          hepatic_vein_T1_C10.l +
          inferior_vena_cava_C12.l +
          inferior_vena_cava_C16.l +
          renal_vein_T1_R18.l +
          inferior_vena_cava_C20.l +
          renal_vein_T1_L22.l +
          inferior_vena_cava_C24.l +
          common_iliac_vein_L56.l +
          common_iliac_vein_R26.l +
          external_iliac_vein_R28.l +
          internal_iliac_vein_T1_R30.l +
          external_iliac_vein_R32.l +
          femoral_vein_R34.l +
          femoral_vein_R38.l +
          profunda_femoris_vein_T2_R40.l +
          femoral_vein_R42.l +
          femoral_vein_R46.l +
          popliteal_vein_R48.l +
          anterior_tibial_vein_T4_R50.l +
          popliteal_vein_R52.l +
          posterior_tibial_vein_T6_R54.l +
          external_iliac_vein_L58.l +
          internal_iliac_vein_T1_L60.l +
          external_iliac_vein_L62.l +
          femoral_vein_L64.l +
          femoral_vein_L68.l +
          profunda_femoris_vein_T2_L70.l +
          femoral_vein_L72.l +
          femoral_vein_L76.l +
          popliteal_vein_L78.l +
          anterior_tibial_vein_T4_L80.l +
          popliteal_vein_L82.l +
          posterior_tibial_vein_T6_L84.l +
          brachiocephalic_vein_R90.l +
          brachiocephalic_vein_L124.l +
          vertebral_vein_R92.l +
          brachiocephalic_vein_R94.l +
          subclavian_vein_R96.l +
          internal_jugular_vein_R122.l +
          external_jugular_vein_R98.l +
          subclavian_vein_R100.l +
          axillary_vein_R102.l +
          brachial_vein_R104.l +
          brachial_vein_R114.l +
          brachial_vein_R108.l +
          ulnar_vein_T7_R110.l +
          brachial_vein_R118.l +
          radial_vein_T3_R120.l +
          vertebral_vein_L126.l +
          brachiocephalic_vein_L128.l +
          subclavian_vein_L130.l +
          internal_jugular_vein_L156.l +
          external_jugular_vein_L132.l +
          subclavian_vein_L134.l +
          axillary_vein_L136.l +
          brachial_vein_L138.l +
          brachial_vein_L148.l +
          brachial_vein_L142.l +
          ulnar_vein_T7_L144.l +
          brachial_vein_L152.l +
          radial_vein_T3_L154.l +
          coronary_veins.l +
          splachnic_vein.l;
      end Systemic_volumes;

      model Systemic_con
        extends ADAN_main.Components.AdanVenousRed.Systemic_volumes;
        replaceable ADAN_main.Components.Auxiliary.AcausalConnector.Pq_terminator_v
          pq_terminator_v(v=-v_aov)
          annotation (Placement(transformation(extent={{-298,78},{-278,98}})));
        replaceable ADAN_main.Components.Auxiliary.AcausalConnector.Pq_terminator_p pq_terminator_sup_vc(u=u_ra)
          annotation (Placement(transformation(extent={{380,160},{360,180}})));
        replaceable ADAN_main.Components.Auxiliary.AcausalConnector.Pq_terminator_p pq_terminator_inf_vc(u=u_ra)
          annotation (Placement(transformation(extent={{480,-10},{460,10}})));
        Real u_ra(unit = "Pa") = port_b.pressure;
        Real v_aov(unit = "m3.s-1") = port_a.q;

        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a port_a annotation (
            Placement(transformation(extent={{-330,70},{-310,90}}),
              iconTransformation(extent={{-330,-10},{-310,10}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b port_b annotation (
            Placement(transformation(extent={{450,70},{470,90}}), iconTransformation(
                extent={{430,-10},{450,10}})));

      equation
        port_a.pressure = pq_terminator_v.u;
        port_b.q = -pq_terminator_inf_vc.v - pq_terminator_sup_vc.v;

        connect(ascending_aorta_A.port_a, pq_terminator_v.port_a) annotation (Line(
            points={{-263,87.5},{-270,87.5},{-270,88},{-278,88}},
            color={0,0,0},
            thickness=1));

        connect(superior_vena_cava_C2.port_b, pq_terminator_sup_vc.port_a)
          annotation (Line(
            points={{317,165.5},{340.5,165.5},{340.5,170},{360,170}},
            color={0,0,0},
            thickness=1));
        connect(inferior_vena_cava_C8.port_b, pq_terminator_inf_vc.port_a)
          annotation (Line(
            points={{429,-2.5},{444.5,-2.5},{444.5,0},{460,0}},
            color={0,0,0},
            thickness=1));
      end Systemic_con;

      model Systemic_baroreflex
          extends ADAN_main.Components.AdanVenousRed.Systemic_con(
            aortic_arch_C46(UseDistentionOutput=true), internal_carotid_R8_A(
              UseDistentionOutput=true));
        ADAN_main.Components.Baroreflex.Baroreflex baroreflex annotation (
            Placement(transformation(extent={{-206,156},{-186,136}})));
        Physiolibrary.Types.RealIO.FractionOutput  phi_baroreflex
                                                      annotation (Placement(
              transformation(extent={{-176,134},{-156,154}}),
                                                          iconTransformation(extent={{-24,164},
                  {-4,184}})));
        Baroreflex.Baroreceptor baroreceptor_aortic(epsilon_start=0.76, s_start=
             0.91) annotation (Placement(transformation(extent={{-236,124},{-216,
                  144}})));
        Baroreflex.Baroreceptor baroreceptor_carotid(epsilon_start=0.4, s_start=
             0.95) annotation (Placement(transformation(extent={{-236,152},{-216,
                  172}})));
      equation
        connect(baroreflex.phi,phi_baroreflex)  annotation (Line(points={{-185.8,146},
                {-174,146},{-174,144},{-166,144}},
                                    color={0,0,127}));
        connect(baroreceptor_carotid.fbr, baroreflex.carotid_BR) annotation (
            Line(points={{-216,162},{-212,162},{-212,156},{-206,156}}, color={0,
                0,127}));
        connect(baroreceptor_aortic.fbr, baroreflex.aortic_BR) annotation (Line(
              points={{-216,134},{-212,134},{-212,136},{-206,136}}, color={0,0,
                127}));
        connect(baroreceptor_aortic.d, aortic_arch_C46.distentionFraction)
          annotation (Line(points={{-236,134},{-242,134},{-242,100},{-126,100},
                {-126,92.5},{-127,92.5}}, color={0,0,127}));
        connect(baroreceptor_carotid.d, internal_carotid_R8_A.distentionFraction)
          annotation (Line(points={{-236,162},{-242,162},{-242,198},{-7,198},{
                -7,194.5}}, color={0,0,127}));
      end Systemic_baroreflex;

      model PulmonaryComponent
        extends ADAN_main.Components.Auxiliary.PulmonaryBase;

        ADAN_main.Components.AdanVenousRed._b580e.Pulmonary pulmonary(
          u_pas(displayUnit="Pa"),
          u_la= port_b.pressure,
          v_puv= port_a.q,
          thoracic_pressure = _thoracic_pressure,
          t=time)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        volume = pulmonary.total_stressed_volume;

          port_a.pressure =pulmonary.u_pas;
          port_b.q =-pulmonary.v_pvn;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                Text(
                extent={{-100,60},{100,100}},
                lineColor={0,0,0},
                textString="ADAN_VR")}),                               Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end PulmonaryComponent;

      model Systemic_simplest
        extends Systemic_interfaces(
          redeclare model Systemic_vein =
              ADAN_main.Components.Vessel_modules.vp_type_tension_based,
            UseThoracic_PressureInput=true,
            UsePhi_Input=true);

        ADAN_main.Components.AdanVenousRed._b580e.Parameters_Venous_cellml.Parameters_Systemic
          Parameters_Systemic1
          annotation (Placement(transformation(extent={{-96,-87},{-76,-82}})));
        replaceable ADAN_main.Components.Vessel_modules.vv_type_thoracic
          ascending_aorta_A constrainedby
          ADAN_main.Components.Vessel_modules.Interfaces.bg_vessel(
          l=Parameters_Systemic1.l_ascending_aorta_A,
          E=Parameters_Systemic1.E_ascending_aorta_A,
          r=Parameters_Systemic1.r_ascending_aorta_A,
          UseInertance=true)
            annotation (Placement(
              transformation(extent={{-263,85},{-243,90}})),
            __Dymola_choicesAllMatching=true);

        Systemic_artery_thoracic ascending_aorta_B(
          UseDistentionOutput=true,
          l=0.0025,
          E = Parameters_Systemic1.E_ascending_aorta_B,
          r=0.007)
        annotation (Placement(transformation(extent={{-238,85},{-218,90}})));

        Systemic_tissue posterior_tibial_T4_R236(
            Ra = tissueParameters.Ra_posterior_tibial_T4_R236,
            Rv = tissueParameters.Rv_posterior_tibial_T4_R236,
            I = tissueParameters.I_posterior_tibial_T4_R236,
            C = tissueParameters.C_posterior_tibial_T4_R236,
            zpv =  tissueParameters.Zpv_posterior_tibial_T4_R236,
          nominal_pressure=tissueParameters.tissue_pressure)
        annotation (Placement(transformation(extent={{35,-5},{55,0}})));
          Systemic_vein inferior_vena_cava_C8(
            l=0.0025,
            E=Parameters_Venous1.E_inferior_vena_cava_C8,
          r=0.007,
          redeclare Vessel_modules.Interfaces.compliance_dataFit1
            compliant_vessel)
            annotation (Placement(transformation(extent={{409,-5},{429,0}})));
        ADAN_main.Components.AdanVenousRed._b580e.Parameters_Venous_cellml.Parameters_Venous
          Parameters_Venous1
          annotation (Placement(transformation(extent={{-69,-87},{-49,-82}})));

        replaceable SystemicTissueParameters.SystemicTissueParameters_Calculated
          tissueParameters(
              total_zpv=0.0027,
              stressed_volume=0.0004,
              qf_posterior_tibial_T4_R236=1) constrainedby
          SystemicTissueParameters.SystemicTissueParameters annotation (
            Placement(transformation(extent={{-140,-98},{-120,-78}})),
            __Dymola_choicesAllMatching=true);
        replaceable ADAN_main.Components.Auxiliary.AcausalConnector.Pq_terminator_v
          pq_terminator_v(v=-v_aov)
          annotation (Placement(transformation(extent={{-298,78},{-278,98}})));
        replaceable ADAN_main.Components.Auxiliary.AcausalConnector.Pq_terminator_p pq_terminator_sup_vc(u=u_ra)
          annotation (Placement(transformation(extent={{380,160},{360,180}})));
        replaceable ADAN_main.Components.Auxiliary.AcausalConnector.Pq_terminator_p pq_terminator_inf_vc(u=u_ra)
          annotation (Placement(transformation(extent={{480,-10},{460,10}})));
        Real u_ra(unit = "Pa") = port_b.pressure;
        Real v_aov(unit = "m3.s-1") = port_a.q;

        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a port_a annotation (
            Placement(transformation(extent={{-330,70},{-310,90}}),
              iconTransformation(extent={{-330,-10},{-310,10}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_b port_b annotation (
            Placement(transformation(extent={{450,70},{470,90}}), iconTransformation(
                extent={{430,-10},{450,10}})));

        Baroreflex.Baroreflex baroreflex annotation (Placement(transformation(
                extent={{-150,152},{-130,132}})));
        Baroreflex.Baroreceptor baroreceptor_aortic(epsilon_start=0.76, s_start=
             0.91) annotation (Placement(transformation(extent={{-180,120},{-160,
                  140}})));
        Baroreflex.Baroreceptor baroreceptor_carotid(epsilon_start=0.4, s_start=
             0.95) annotation (Placement(transformation(extent={{-180,148},{-160,
                  168}})));
        Physiolibrary.Types.RealIO.FractionOutput  phi_baroreflex
                                                      annotation (Placement(
              transformation(extent={{-120,130},{-100,150}}),
                                                          iconTransformation(extent={{-24,164},
                  {-4,184}})));
      equation
        port_a.pressure = pq_terminator_v.u;
        port_b.q = -pq_terminator_inf_vc.v - pq_terminator_sup_vc.v;

        connect(ascending_aorta_A.port_a, pq_terminator_v.port_a) annotation (Line(
            points={{-263,87.5},{-270,87.5},{-270,88},{-278,88}},
            color={0,0,0},
            thickness=1));

        connect(inferior_vena_cava_C8.port_b, pq_terminator_inf_vc.port_a)
          annotation (Line(
            points={{429,-2.5},{444.5,-2.5},{444.5,0},{460,0}},
            color={0,0,0},
            thickness=1));


        connect(ascending_aorta_B.port_a,ascending_aorta_A.port_b) annotation (Line(points={{-238,
                87.5},{-243,87.5}},                                                                                        thickness=1));
        connect(ascending_aorta_B.port_b, posterior_tibial_T4_R236.port_a)
          annotation (Line(
            points={{-218,87.5},{-90,87.5},{-90,-2.5},{35,-2.5}},
            color={0,0,0},
            thickness=1));
        connect(posterior_tibial_T4_R236.port_b, inferior_vena_cava_C8.port_a)
          annotation (Line(
            points={{55,-2.5},{232.5,-2.5},{232.5,-2.5},{409,-2.5}},
            color={0,0,0},
            thickness=1));
        connect(baroreflex.phi,phi_baroreflex)  annotation (Line(points={{-129.8,142},
                {-118,142},{-118,140},{-110,140}},
                                    color={0,0,127}));
        connect(baroreceptor_carotid.fbr,baroreflex. carotid_BR) annotation (
            Line(points={{-160,158},{-156,158},{-156,152},{-150,152}}, color={0,
                0,127}));
        connect(baroreceptor_aortic.fbr,baroreflex. aortic_BR) annotation (Line(
              points={{-160,130},{-156,130},{-156,132},{-150,132}}, color={0,0,
                127}));
        connect(baroreceptor_aortic.d, ascending_aorta_B.distentionFraction)
          annotation (Line(points={{-180,130},{-228,130},{-228,92.5}}, color={0,0,127}));
        connect(baroreceptor_carotid.d, ascending_aorta_B.distentionFraction)
          annotation (Line(points={{-180,158},{-228,158},{-228,92.5}}, color={0,0,127}));
          annotation (Diagram(coordinateSystem(extent={{-320,-100},{440,200}})),
              Icon(coordinateSystem(extent={{-320,-100},{440,200}}), graphics={
                        Bitmap(extent={{-320,-100},{440,200}}, fileName=
                    "modelica://Physiolibrary/Resources/Icons/perfusion.png")}));
      end Systemic_simplest;
    end AdanVenousRed;

    model MyDelay
      extends Modelica.Blocks.Interfaces.SISO;
      parameter Modelica.SIunits.Time delayTime(start=1)
        "Delay time of output with respect to input signal";

    //  parameter Physiolibrary.Types.Time phiDelay = 0;
    //  Physiolibrary.Types.RealIO.FractionInput phi_in annotation (Placement(
    //       transformation(extent={{-100,-100},{-60,-60}}),
    //                                                     iconTransformation(extent={{-100,60},
    //               {-60,100}})));

    equation
        if delayTime > 0 then
          y = delay(u, delayTime,  10);
        else
          y = u;
        end if;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Text(
          extent={{8.0,-142.0},{8.0,-102.0}},
          textString="delayTime=%delayTime"),
                                    Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
        Line(
          points={{-92.0,0.0},{-80.7,34.2},{-73.5,53.1},{-67.1,66.4},{-61.4,74.6},{-55.8,79.1},{-50.2,79.8},{-44.6,76.6},{-38.9,69.7},{-33.3,59.4},{-26.9,44.1},{-18.83,21.2},{-1.9,-30.8},{5.3,-50.2},{11.7,-64.2},{17.3,-73.1},{23.0,-78.4},{28.6,-80.0},{34.2,-77.6},{39.9,-71.5},{45.5,-61.9},{51.9,-47.2},{60.0,-24.8},{68.0,0.0}},
          color={0,0,127},
          smooth=Smooth.Bezier),
        Line(
          points={{-62.0,0.0},{-50.7,34.2},{-43.5,53.1},{-37.1,66.4},{-31.4,74.6},{-25.8,79.1},{-20.2,79.8},{-14.6,76.6},{-8.9,69.7},{-3.3,59.4},{3.1,44.1},{11.17,21.2},{28.1,-30.8},{35.3,-50.2},{41.7,-64.2},{47.3,-73.1},{53.0,-78.4},{58.6,-80.0},{64.2,-77.6},{69.9,-71.5},{75.5,-61.9},{81.9,-47.2},{90.0,-24.8},{98.0,0.0}},
          color={160,160,164},
          smooth=Smooth.Bezier)}),                                   Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end MyDelay;

    package Smith

      model PulmonarySmith
        extends Auxiliary.PulmonaryBase(UseThoracic_PressureInput=true);
          Physiolibrary.Hydraulic.Components.ElasticVesselElastance
                                 pulmonaryArteries(
          ZeroPressureVolume=0,
          useExternalPressureInput=true,
          volume_start=3.904e-05,
          Elastance(displayUnit="Pa/m3") = 49195960.956135)
          annotation (Placement(transformation(extent={{-44,-10},{-24,10}})));
        Physiolibrary.Hydraulic.Components.Resistor
                 Rpul(Resistance(displayUnit="(mmHg.s)/ml") = 20691634.526808)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={0,0})));
      Physiolibrary.Hydraulic.Components.ElasticVesselElastance
                             pulmonaryVeins(
          ZeroPressureVolume=0,
          useExternalPressureInput=true,
          volume_start=0.0008269,
          Elastance(displayUnit="Pa/m3") = 973253.4281295)
          annotation (Placement(transformation(extent={{32,-10},{52,10}})));
        Physiolibrary.Types.Constants.PressureConst pressure(k=0) if not
          UseThoracic_PressureInput
          annotation (Placement(transformation(extent={{-32,28},{-24,36}})));
      equation
        volume = pulmonaryArteries.volume + pulmonaryVeins.volume;

        connect(pulmonaryArteries.q_in,Rpul. q_in) annotation (Line(
            points={{-34,2.22045e-16},{-10,2.22045e-16},{-10,0}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None));
        connect(Rpul.q_out,pulmonaryVeins. q_in) annotation (Line(
            points={{10,2.22045e-16},{10,0},{42,0},{42,2.22045e-16}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None));
        connect(thoracic_pressure, pulmonaryVeins.externalPressure)
          annotation (Line(points={{-8,-100},{50,-100},{50,8}}, color={0,0,127}));
        connect(thoracic_pressure, pulmonaryArteries.externalPressure) annotation (
            Line(points={{-8,-100},{50,-100},{50,8},{-26,8}}, color={0,0,127}));
        connect(pulmonaryArteries.q_in, port_a) annotation (Line(
            points={{-34,2.22045e-16},{-36,2.22045e-16},{-36,0},{-100,0}},
            color={0,0,0},
            thickness=1));
        connect(pulmonaryVeins.q_in, port_b) annotation (Line(
            points={{42,0},{100,0}},
            color={0,0,0},
            thickness=1));
        connect(pressure.y, pulmonaryArteries.externalPressure) annotation (Line(
              points={{-23,32},{-20,32},{-20,8},{-26,8}}, color={0,0,127}));
        connect(pressure.y, pulmonaryVeins.externalPressure) annotation (Line(
              points={{-23,32},{18,32},{18,8},{50,8}}, color={0,0,127}));
        annotation (Icon(graphics={Text(
                extent={{-100,60},{100,100}},
                lineColor={0,0,0},
                textString="Smith")}));
      end PulmonarySmith;

      model HeartSmith
        extends Auxiliary.HeartBase;
        parameter Boolean UsePhiInput = false;

        IdealValveResistanceWithMeasurements
                             aorticValve(Pknee=0, _Ron(displayUnit="(mmHg.s)/ml")=
               2399802.97347)
          annotation (Placement(transformation(extent={{-68,-30},{-88,-10}})));
        Physiolibrary.Hydraulic.Components.IdealValveResistance
                             tricuspidValve(Pknee=0, _Ron(displayUnit=
                "(mmHg.s)/ml") = 3159740.5817355)
          annotation (Placement(transformation(extent={{-62,24},{-42,44}})));
        Physiolibrary.Hydraulic.Components.Inertia
                Lav(I(displayUnit="mmHg.s2/ml") = 16250.665802014,
            volumeFlow_start(displayUnit="m3/s") = -1.4e-8) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-44,-20})));
        Physiolibrary.Hydraulic.Components.Inertia
                Lpv(I(displayUnit="mmHg.s2/ml") = 19822.372560862,
            volumeFlow_start(displayUnit="m3/s") = -1.9e-9)
          annotation (Placement(transformation(extent={{32,24},{52,44}})));
        Physiolibrary.Hydraulic.Components.IdealValveResistance
                             pulmonaryValve(Pknee=0, _Ron(displayUnit=
                "(mmHg.s)/ml") = 733273.1307825)
          annotation (Placement(transformation(extent={{62,24},{82,44}})));
        Physiolibrary.Hydraulic.Components.IdealValveResistance
                             mitralValve(Pknee=0, _Ron(displayUnit="(mmHg.s)/ml")=
               2106493.721157)
          annotation (Placement(transformation(extent={{52,-30},{32,-10}})));
        Physiolibrary.Hydraulic.Components.Inertia
                Ltc(I(displayUnit="mmHg.s2/ml") = 10678.18997523,
            volumeFlow_start(displayUnit="m3/s") = 0.0001372)
          annotation (Placement(transformation(extent={{-88,24},{-68,44}})));
        Physiolibrary.Hydraulic.Components.Inertia
                Lmt(I(displayUnit="mmHg.s2/ml") = 10261.557514558,
            volumeFlow_start(displayUnit="m3/s") = 0.0001141) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={70,-20})));
        Smith_VentricularInteraction_flat smith_VentricularInteraction_flat(alphaE=
              settings.heart_alphaE)
          annotation (Placement(transformation(extent={{-14,-10},{16,28}})));
        Physiolibrary.Types.RealIO.FractionInput phi if UsePhiInput
          annotation (Placement(transformation(extent={{-120,50},{-80,90}}),
              iconTransformation(extent={{-120,40},{-80,80}})));
        Physiolibrary.Types.Constants.FractionConst phi0(k=settings.phi0) if not UsePhiInput
          annotation (Placement(transformation(extent={{-88,82},{-80,90}})));
        outer Settings settings
          annotation (Placement(transformation(extent={{40,80},{60,100}})));
      equation
        connect(Lav.q_out,aorticValve. q_in) annotation (Line(
            points={{-54,-20},{-68,-20}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None));
        connect(Ltc.q_out,tricuspidValve. q_in) annotation (Line(
            points={{-68,34},{-62,34}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None));
        connect(Lpv.q_out,pulmonaryValve. q_in) annotation (Line(
            points={{52,34},{62,34}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None));
        connect(mitralValve.q_in,Lmt. q_out) annotation (Line(
            points={{52,-20},{60,-20}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None));
        connect(mitralValve.q_out,Lav. q_in) annotation (Line(
            points={{32,-20},{-34,-20}},
            color={0,0,0},
            thickness=1,
            smooth=Smooth.None));
        connect(sa, aorticValve.q_out) annotation (Line(
            points={{100,100},{100,-60},{-96,-60},{-96,-20},{-88,-20}},
            color={0,0,0},
            thickness=1));
        connect(Lmt.q_in, pv) annotation (Line(
            points={{80,-20},{80,-82},{-100,-82},{-100,-100}},
            color={0,0,0},
            thickness=1));
        connect(sv, Ltc.q_in) annotation (Line(
            points={{-100,100},{-100,34},{-88,34}},
            color={0,0,0},
            thickness=1));
        connect(pulmonaryValve.q_out, pa) annotation (Line(
            points={{82,34},{92,34},{92,-100},{100,-100}},
            color={0,0,0},
            thickness=1));
        connect(tricuspidValve.q_out, smith_VentricularInteraction_flat.rvflow)
          annotation (Line(
            points={{-42,34},{-20,34},{-20,28},{0.7,28}},
            color={0,0,0},
            thickness=1));
        connect(Lpv.q_in, smith_VentricularInteraction_flat.rvflow) annotation (Line(
            points={{32,34},{16,34},{16,28},{0.7,28}},
            color={0,0,0},
            thickness=1));
        connect(Lav.q_in, smith_VentricularInteraction_flat.lvflow) annotation (Line(
            points={{-34,-20},{2,-20},{2,-10},{1,-10}},
            color={0,0,0},
            thickness=1));
        connect(mitralValve.q_out, smith_VentricularInteraction_flat.lvflow)
          annotation (Line(
            points={{32,-20},{2,-20},{2,-10},{1,-10}},
            color={0,0,0},
            thickness=1));
        connect(smith_VentricularInteraction_flat.HR, frequency_input) annotation (
            Line(points={{-11,9},{-55.5,9},{-55.5,0},{-106,0}}, color={0,0,127}));
        connect(smith_VentricularInteraction_flat.Pth, thoracic_pressure_input)
          annotation (Line(points={{13.3,9},{13.3,-44.5},{-8,-44.5},{-8,-100}}, color=
               {0,0,127}));
        connect(phi, smith_VentricularInteraction_flat.phi) annotation (Line(points={{
                -100,70},{-56,70},{-56,24.2},{-11,24.2}}, color={0,0,127}));
        connect(phi0.y, smith_VentricularInteraction_flat.phi) annotation (Line(
              points={{-79,86},{-20,86},{-20,24.2},{-11,24.2}}, color={0,0,127}));
        connect(smith_VentricularInteraction_flat.Pth, P0.y) annotation (Line(
              points={{13.3,9},{13.3,-44.5},{-8,-44.5},{-8,-74},{-4,-74},{-4,
                -75},{9.4369e-16,-75}}, color={0,0,127}));
        connect(smith_VentricularInteraction_flat.HR, HR0.y) annotation (Line(
              points={{-11,9},{-55.5,9},{-55.5,0},{-75,0}}, color={0,0,127}));
        annotation (Diagram(graphics={  Rectangle(extent={{-84,54},{92,-46}},
                lineColor={28,108,200})}), Icon(graphics={Text(
                extent={{-100,60},{100,100}},
                lineColor={0,0,0},
                textString="Smith")}));
      end HeartSmith;

      model Smith_VentricularInteraction_flat
          import Physiolibrary.Types.*;
        Volume Vsept(start=V0sept), Vrv(start=0.0001042), Vlv(start=
              0.0001042), Vperi;
        parameter Volume V0sept=2e-06, V0peri = 0.0002;

        Pressure Psept, Pperi;
        parameter Pressure Pi0sept=148.00118226939, Pi0rv=28.757638965416, Pi0lv=16.038683206025, Pi0peri=66.701190423724
          "peak isovolumic pressure";

        parameter HydraulicElastance Essept0 = 6499999676.0309, Esrv0= 77993596.637775, Eslv0=383941811.27772;

      //   parameter Physiolibrary.Types.Volume VS0sept = 1e-3 "Volume Threshold for linear Frank-starling effect";
      //   parameter Physiolibrary.Types.Volume VS0lv = 1e-3 "Volume Threshold for linear Frank-starling effect";
      //   parameter Physiolibrary.Types.Volume VS0rv = 1e-3 "Volume Threshold for linear Frank-starling effect";

      //   parameter Real Escale = 1;
        Physiolibrary.Types.Fraction PhiEffect = (1 + alphaE*(phi - phi0));

        HydraulicElastance Essept = Essept0*PhiEffect;
        HydraulicElastance Esrv = Esrv0*PhiEffect;
        HydraulicElastance Eslv = Eslv0*PhiEffect    "elastance of systole";
        // HydraulicElastance Essept = Essept0*(1 + alphaE*(phi - phi0))*(1 + (Escale*Vsept/VS0sept-1)*(tanh(Vsept/VS0sept-2)+1)/2);
        // HydraulicElastance Esrv = Esrv0*(1 + alphaE*(phi - phi0))*(1 + (Escale*Vrv/VS0rv-1)*(tanh(Vrv/VS0rv-2)+1)/2);
        // HydraulicElastance Eslv = Eslv0*(1 + alphaE*(phi - phi0))*(1 + (Escale*Vlv/VS0lv-1)*(tanh(Vlv/VS0lv-2)+1)/2)
        //   "elastance of systole";
        parameter Real A=1 "Multiplier of driving function";
        Real B= 60/HP;
        // Real CC=HP*SystolicFraction;
        Real CC=HP*ts;
      //  parameter Physiolibrary.Types.Fraction SystolicFraction = 0.5;
          Time tm;
          discrete Time HP "heart period";
          discrete Time t0 "time of beginning of the cardiac cycle";
          discrete Time ts "duration of systole";
          Time td = HP - ts "duration of diastole";
          parameter Real lambdas= 435000 "Lambda of septum [1/m3]";
          parameter Real lambdarv = 23000 "Lambda of RV [1/m3]";
          parameter Real lambdalv = 33000 "Lambda of LV [1/m3]";
          parameter Real lambdaperi = 30000 "Lambda of pericardium [1/m3]";
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a rvflow annotation (
           Placement(transformation(extent={{-48,20},{-28,40}}),
              iconTransformation(extent={{-12,90},{8,110}})));
        Physiolibrary.Hydraulic.Interfaces.HydraulicPort_a lvflow annotation (
           Placement(transformation(extent={{-46,-22},{-26,-2}}),
              iconTransformation(extent={{-10,-110},{10,-90}})));
        RealIO.FrequencyInput HR annotation (Placement(transformation(extent=
                  {{-78,-40},{-38,0}}), iconTransformation(extent={{-100,-20},
                  {-60,20}})));
        RealIO.PressureInput Pth annotation (Placement(transformation(extent=
                  {{-6,24},{14,44}}), iconTransformation(
              extent={{-20,-20},{20,20}},
              rotation=180,
              origin={82,0})));
        Physiolibrary.Types.RealIO.FractionInput
                              phi
                                 annotation (Placement(transformation(extent={{-80,-100},
                  {-40,-60}}),          iconTransformation(extent={{-100,60},{-60,100}})));
        Physiolibrary.Types.Fraction phi0 = 0.25;
        parameter Physiolibrary.Types.Fraction alphaE = 0 "linear dependency of elastances (Essept, Esrv, Eslv) on phi";
        Physiolibrary.Types.Pressure Prv = rvflow.pressure - Pperi;
        Physiolibrary.Types.Pressure Plv = lvflow.pressure - Pperi;
        parameter Physiolibrary.Types.Fraction alphaDriving = 0 "Experimental sensitivity of driving funtion on phi";
        Real driving = (1 + alphaDriving*(phi - phi0))*A*exp(-B*(tm - CC)^2) "Linear dependency of driving on phi";

        Modelica.Blocks.Interfaces.BooleanOutput
                              beat = tm >= pre(HP)
                                 annotation (Placement(transformation(extent={{-28,-120},
                  {12,-80}}),           iconTransformation(extent={{-60,-80},{-100,-40}})));
                  Physiolibrary.Types.Volume ESV_LV;
                  Physiolibrary.Types.Volume ESV_RV;
                  Physiolibrary.Types.Volume EDV_LV;
                  Physiolibrary.Types.Volume EDV_RV;
                  Physiolibrary.Types.Volume SV_LV;
                  Physiolibrary.Types.Volume SV_RV;
                  Physiolibrary.Types.VolumeFlowRate CO_LV = SV_LV/HP;
                  Physiolibrary.Types.VolumeFlowRate CO_RV = SV_RV/HP;
                  parameter Time ts_a1 = 0.1 "Part of calculation for t_systole = a1 + a2*period. Guessed from Bombardino 2008 doi:10.1186/1476-7120-6-15";
                  parameter Real ts_a2 = 0.2 "Part of calculation for t_systole = a1 + a2*period. Guessed from Bombardino 2008 doi:10.1186/1476-7120-6-15";
      equation
        //timing
        tm = time - pre(t0);
        when {initial(),beat} then
          HP = 1/HR;
          t0 = time;
          //    ts = 0.16 + 0.3*HP;
          ts = ts_a1 + ts_a2*HP "Guessed from Bombardino 2008 doi:10.1186/1476-7120-6-15";
          EDV_LV = Vlv;
          EDV_RV = Vrv;
        end when;

        when time > t0 + ts then
          ESV_LV = Vlv;
          ESV_RV = Vrv;
          SV_LV = EDV_LV - ESV_LV;
          SV_RV = EDV_RV - ESV_RV;
        end when;
        //  septum
        Psept = lvflow.pressure - rvflow.pressure;
        Psept = (Vsept - V0sept)*driving*Essept +
          (1 - driving)*Pi0sept*(exp(lambdas*Vsept) - 1);
        // rightventricle
        rvflow.pressure - Pperi = (Vrv + Vsept)*driving*Esrv +
          (1 - driving)*Pi0rv*(exp(lambdarv*(Vrv + Vsept)) - 1);
        der(Vrv) = rvflow.q;
        //leftventricle
        lvflow.pressure - Pperi = (Vlv - Vsept)*driving*Eslv +
          (1 - driving)*Pi0lv*(exp(lambdalv*(Vlv - Vsept)) - 1);
        der(Vlv) = lvflow.q;
        //pericardium
        Vperi = Vrv + Vlv;
        Pperi = Pth + Pi0peri*(exp(lambdaperi*(Vperi - V0peri)) - 1);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={
                  {-100,-100},{100,100}}), graphics={Text(
                extent={{102,-32},{76,-20}},
                lineColor={0,0,255},
                fillColor={255,170,170},
                fillPattern=FillPattern.Forward,
                        textString="Pth"),Text(
                extent={{-100,-22},{-74,-34}},
                lineColor={0,0,255},
                        textString="HR"),Rectangle(
                        extent={{-20,80},{20,-60}},
                        lineColor={0,0,255},
                fillPattern=FillPattern.Solid,
                        fillColor={0,0,255}),Text(
                extent={{-100,-60},{100,-80}},
                lineColor={0,0,255},
                textString="%name"),      Text(
                extent={{-64,70},{-38,58}},
                lineColor={0,0,255},
                textString="phi")}));
      end Smith_VentricularInteraction_flat;

      model IdealValveResistanceWithMeasurements
        "Adds additional output calculations"
        extends Physiolibrary.Hydraulic.Components.IdealValveResistance;

        Physiolibrary.Types.Pressure BPp;
        Physiolibrary.Types.Pressure BPs;
        Physiolibrary.Types.Pressure BPd;


        Physiolibrary.Types.Pressure BP_max(start = 0);
        Physiolibrary.Types.Pressure BP_min(start = 200*133);

        Physiolibrary.Types.Pressure BPao = q_out.pressure;

        parameter Physiolibrary.Types.Time tau = 1e-3;
      equation
      //   if open then
      //     BP_max = max(BPao, BP_max);
      //   else
      //     BP_max = 0;
      //   end if;

        if open and BP_max < BPao then
          der(BP_max)*tau = BPao - BP_max;
        else
          der(BP_max) = 0;
        end if;

        if not open and BP_min > BPao then
          der(BP_min)*tau = BPao - BP_min;
        else
          der(BP_min) = 0;
        end if;

      // workaround for OM
        // when open then
      when passableVariable > Modelica.Constants.eps then
          BPd = pre(BP_min);
          reinit(BP_min, 200*133);
          BPp = BPs - BPd;
        end when;

        // workaround for OM
      //  when not open then
        when passableVariable < - Modelica.Constants.eps then
          BPs = pre(BP_max);
          reinit(BP_max, 0);
      //     BPp = BPs - BPd;
        end when;


      end IdealValveResistanceWithMeasurements;
    end Smith;

    package Baroreflex

      model Baroreceptor
        Physiolibrary.Types.RealIO.FractionInput d "The distension ratio r/r0. Should be around 1, but not necesarily exactly 1, as it is compensated by other paraemters"
        annotation (Placement(transformation(
                extent={{-110,-8},{-90,12}}),  iconTransformation(extent={{-120,
                  -20},{-80,20}})));
        Modelica.Blocks.Interfaces.RealOutput fbr( unit = "Hz") = f0*s*(delta/(delta + delta0)) "Baroreceptor firing frequency" annotation (Placement(transformation(
                extent={{84,-10},{104,10}}),   iconTransformation(extent={{80,-20},
                  {120,20}})));

        Real epsilon( start = epsilon_start) "Averaged distension ratio";
        parameter Physiolibrary.Types.Time Ts = 30 "Time constant for averaging";
        Real delta=max(d - epsilon, 0) "Positive peaks detected";
        parameter Real f0( unit = "Hz")= 300 "Base firing frequency";
        parameter Real delta0 = 0.4965 "Baseline delta";

        Real s(start = s_start);
        parameter Real a(unit="s-1") = 0.0651;
        parameter Real b(unit="s-1") = 0.2004;
        parameter Real epsilon_start = 1.075;
        parameter Real s_start = 0.85;
        parameter Modelica.SIunits.Time resetAt = 0 "resets initial conditions to counter transients";

      equation

        when time > resetAt then
          reinit(epsilon, epsilon_start);
          reinit(s, s_start);
        end when;

        der(epsilon) =(d - epsilon)/Ts;
        der(s) =a*(1 - s) - b*s*(delta/(delta + delta0));

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Baroreceptor;

      model Baroreflex
        Modelica.Blocks.Interfaces.RealInput carotid_BR annotation (Placement(transformation(
                extent={{-118,-68},{-78,-28}}),
                                              iconTransformation(extent={{-120,-120},{
                  -80,-80}})));
        Modelica.Blocks.Interfaces.RealInput aortic_BR annotation (Placement(transformation(
                extent={{-114,48},{-74,88}}), iconTransformation(extent={{-120,80},{-80,
                  120}})));

      Real fiSN(start = fiSN_start);
      parameter Real fsn( unit = "s-1") = 0.041;
      parameter Real f1 = 0.0046;
      parameter Real g = 0.66;
      Real aorticWeight = 2*g*aortic_BR;
      Real carotidWeight = 2*(1-g)*carotid_BR;
      parameter Modelica.SIunits.Time resetAt = -1;
      parameter Real fiSN_start = 0.25;
        Physiolibrary.Types.RealIO.FractionOutput phi = fiSN
                                                      annotation (Placement(
              transformation(extent={{96,70},{116,90}}),  iconTransformation(extent={{92,-10},
                  {112,10}})));
      equation
        when time > resetAt then
          reinit(fiSN, fiSN_start);
        end when;

        der(fiSN) = fsn*(1-fiSN) - fiSN*f1*(aorticWeight + carotidWeight);
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Baroreflex;

      model HeartRate

        parameter Physiolibrary.Types.Frequency H0=0.46666666666667;
        parameter Physiolibrary.Types.Frequency H1=2.6;

        Physiolibrary.Types.RealIO.FrequencyOutput HR annotation (Placement(
              transformation(extent={{90,-10},{110,10}}), iconTransformation(extent={{92,-10},
                  {112,10}})));
        Physiolibrary.Types.RealIO.FractionInput
                                             phi annotation (Placement(transformation(
                extent={{-120,-20},{-80,20}}),iconTransformation(extent={{-120,-20},{-80,
                  20}})));

      equation
        HR = H0 + H1*phi;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end HeartRate;
    end Baroreflex;

    model Settings
      import Physiolibrary.Types.*;

      // general
      parameter Fraction phi0=0.25   "Baseline resting phi";

      // HEART
      parameter Physiolibrary.Types.Fraction heart_alphaE = 0 "linear dependency of active elastance on phi" annotation(Dialog(group = "Heart"));

      // arteries
      parameter Pressure tissues_nominal_arterioles_pressure=13332.2387415
                                                                   "Nominal arterial pressure for initialization and calculation of arterioles resistance" annotation(Dialog(group = "Arteries"));
      parameter Boolean arteries_UseVasoconstrictionEffect=false "Change compliance of large arteries by phi and exercise" annotation(choices(checkBox=true), Dialog(group = "Arteries"));
      parameter Fraction exercise_factor_on_arterial_compliance = 0 "Effect of venoconstriction"  annotation(Dialog(group = "Arteries", enable = arteries_UseVasoconstrictionEffect));

      parameter Fraction R_vc = 0 "Effect fraction of venoconstriction on compliance" annotation(Dialog(group = "Arteries", enable = arteries_UseVasoconstrictionEffect));

    // tissues
      parameter Pressure tissues_nominal_pressure=2666.4477483 "Tissues and small arteries and veins nominal pressure. Used for calculation of arteriole and venule resistances and tissues compliance"
                                                       annotation(Dialog(group = "Tissues"));
      parameter Volume tissues_nominal_zpv=0.002304 "Tissues and small arteries and veins nominal pressure. Used for calculation of tissues total volume" annotation(Dialog(group = "Tissues"));
      parameter Volume tissues_nominal_stressed_volume=0.000795 "Tissues and small arteries and veins nominal stressed volume. Used for calculation of nominal tissues compliance"  annotation(Dialog(group = "Tissues"));
      parameter VolumeFlowRate tissues_nominal_cardiac_output=9.98e-05 "Nominal flow through systemic tissues. Used for calculation of arteriole and venule resistances."
                                                                      annotation(Dialog(group = "Tissues"));

      parameter Fraction Ra_factor=1   "Exponential factor affecting arterioles resistance" annotation(Dialog(group = "Tissues"));
      parameter Time tissues_Ra_tau=1.2e-05 "guess from Pstras, Math Med Biol 2017" annotation(Dialog(group = "Tissues"));
      parameter Fraction Rv_factor = 0 "Exponential factor affecting venules resistance" annotation(Dialog(group = "Tissues"));

      parameter Boolean UseNonLinear_TissuesCompliance = false "Use nonlinear tissues PV chars, using the log Vmax suggested by Hardy and collins and used in Pstras 2017" annotation(choices(checkBox=true), Dialog(group = "Tissues"));
      parameter Boolean tissues_UseStraighteningReaction2Phi = false "Use V_n dependency on phi (instead of V_u and V_max as used in Pstras)" annotation(choices(checkBox=true), Dialog(group = "Tissues"));
    //  parameter Boolean UseNonLinear_TissuesCompliance_PhiEffect = false annotation(choices(checkBox=true), Dialog(enable = UseNonLinear_TissuesCompliance, group = "Tissues"));
      parameter Fraction tissuesCompliance_PhiEffect = 0 "Effect on tissue's compliance" annotation(Dialog(enable = UseNonLinear_TissuesCompliance, group = "Tissues"));
      parameter Fraction exercise_factor_on_tissue_compliance = 0 "Effect of venoconstriction"  annotation(Dialog(group = "Tissues", enable = UseNonLinear_TissuesCompliance));

      parameter Fraction tissues_gamma=1   "Nonlinear tissues compliance steepness to set Vmax. Vmax = Vn + gamma*(Vn - zpv)"
        annotation (Dialog(enable = UseNonLinear_TissuesCompliance, group = "Tissues"));
      parameter Fraction exercise_factor = 1 annotation (Dialog(group = "Tissues"));
        //  parameter Fraction tissues_compliance_phi_shift=0 "complinace dependency on phi" annotation(Dialog(group = "Tissues"));


    // veins
      parameter Pressure tissues_nominal_venules_pressure=666.611937075
                                                                "Venous pressure used for initialization and calculation of tissue venules resistances"
        annotation(Dialog(group = "Veins"));
      parameter Boolean UseNonLinear_VenousCompliance = false "Use nonlinear compliance model in all large veins" annotation(choices(checkBox=true), Dialog(group = "Veins"));
      parameter Boolean veins_UsePhiEffect=true "Use venoconstriction in all large veins" annotation(choices(checkBox=true), Dialog(group = "Veins"));
    //  parameter Pressure venous_p0 = 666.6 "Venous pressure used for initialization"    annotation(Dialog(group = "Veins"));
      // parameter Real gamma;
      // parameter Real alpha;
      parameter Fraction veins_gamma=0.5    "Fraction of minimal collapsing diameter to nominal diameter"    annotation(Dialog(group = "Veins"));
      parameter Fraction veins_alpha=5   "ONLY for unused alpha-based tension model> how many times the tension is larger for maximal activation from resting activation at nominal diameter"    annotation(Dialog(group = "Veins"));
      parameter Time veins_activation_tau(displayUnit="s")=0.1
                                                "Integration delay for venous tone activation"    annotation(Dialog(group = "Veins"));
      parameter Fraction venous_diameter_correction=1.5 "Venous diameter enlargement to correspond to the arterial side" annotation(Dialog(group = "Veins"));

      parameter Boolean veins_ignoreViscosityResistance = true "Ignore viscoelastic resistance in all large veins (negligible, but computation-heavy)" annotation(choices(checkBox=true), Dialog(group = "Veins"));
      parameter Boolean veins_limitExternalPressure =  true "limit Pext in large veins, when they are already empty" annotation(choices(checkBox=true), Dialog(group = "Veins"));
    //  parameter Boolean veins_UseViscoElasticDelay = false;
    //   parameter Fraction gamma =   0.5

    //

      annotation (defaultComponentName =     "settings",
                 defaultComponentPrefixes = "inner",
                 missingInnerMessage =      "The Settings object is missing. Add a Settings to the top level and make it 'outer'",
      Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
              points={{80,100},{40,100},{40,98},{40,60},{-60,-40},{-100,-40},{-100,-80},
                  {-80,-60},{-60,-80},{-80,-100},{-40,-100},{-40,-60},{60,40},{100,40},
                  {100,80},{80,60},{60,80},{80,100}},
              lineColor={0,0,0},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid), Polygon(
              points={{-80,100},{-40,100},{-40,98},{-40,60},{60,-40},{100,-40},{100,
                  -80},{80,-60},{60,-80},{80,-100},{40,-100},{40,-60},{-60,40},{-100,
                  40},{-100,80},{-80,60},{-60,80},{-80,100}},
              lineColor={0,0,0},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid)}),                      Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Settings;
  end Components;

  package AdanVenousRed_Safaei

    package Experiments

      model simple_phi_base
        extends CVS_7af(
          redeclare Components.Smith.PulmonarySmith pulmonaryComponent(
              UseThoracic_PressureInput=true),
          redeclare Components.Smith.HeartSmith heartComponent(UseFrequencyInput=
                true, UseThoracicPressureInput=true),
          redeclare Components.AdanVenousRed.Systemic_simplest Systemic1(
            UseThoracic_PressureInput=true,
            UsePhi_Input=true,
            baroreflex(resetAt=-1,
              fsn=0.025,
              g=1.0),
            baroreceptor_aortic(epsilon_start=1.19, delta0=0.3),
            baroreceptor_carotid(epsilon_start=1.06, s_start=0.96),
            ascending_aorta_B(
              UseInertance=false,
              UseVasoconstrictionEffect=false,
              l(displayUnit="m") = 2.5),
            inferior_vena_cava_C8(
              l=2.5,
              r(displayUnit="mm") = 0.008,
              compliant_vessel(tau(displayUnit="d") = 604800))),
          phi(
            amplitude=0.74,
            rising=200,
            width=50,
            period=1000,
            startTime=50));
      equation
        connect(phi.y, condPhi.u) annotation (Line(points={{-49,4},{-38,4},{-38,
                4.00002},{-27.2,4.00002}}, color={0,0,127}));
        connect(phi.y, condHR.u) annotation (Line(points={{-49,4},{-38,4},{-38,
                60},{40.8,60}}, color={0,0,127}));
        annotation (experiment(
            StopTime=400,
            Interval=0.11,
            Tolerance=1e-05,
            __Dymola_Algorithm="Cvode"));
      end simple_phi_base;
    end Experiments;

    partial model CVS_7af
      replaceable Components.AdanVenousRed.Systemic_baroreflex Systemic1(
                baroreceptor_aortic(delta0=0.6,
          epsilon_start=1.23,
          s_start=0.922),
                baroreceptor_carotid(delta0=0.3,
          epsilon_start=1.07,
          s_start=0.945),
                baroreflex(fsn=0.021),
        UseThoracic_PressureInput=true,
        UsePhi_Input=true,
        redeclare model Systemic_vein =
            ADAN_main.Components.Vessel_modules.vp_type_tension_based,
        redeclare
          Components.AdanVenousRed.SystemicTissueParameters.SystemicTissueParameters_Calculated
          tissueParameters,
        femoral_vein_R34(LimitBackflow=true),
        femoral_vein_L64(LimitBackflow=true),
        superior_vena_cava_C88(UseOuter_thoracic_pressure=true),
        superior_vena_cava_C2(UseOuter_thoracic_pressure=true),
        inferior_vena_cava_C8(UseOuter_thoracic_pressure=true),
        hepatic_vein_T1_C10(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        splachnic_vein(UseOuter_thoracic_pressure=true, thoracic_pressure_ratio=
             0.8),
        renal_vein_T1_R18(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        internal_iliac_vein_T1_R30(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        internal_iliac_vein_T1_L60(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        renal_vein_T1_L22(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        inferior_vena_cava_C20(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        inferior_vena_cava_C16(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        inferior_vena_cava_C12(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        abdominal_aorta_C114(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        abdominal_aorta_C136(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        abdominal_aorta_C164(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        abdominal_aorta_C176(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        abdominal_aorta_C188(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        abdominal_aorta_C192(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        mesenteric_artery(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        common_iliac_R216(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8),
        common_iliac_L194(UseOuter_thoracic_pressure=true,
            thoracic_pressure_ratio=0.8)) constrainedby
        Components.AdanVenousRed.Systemic_interfaces
                                               annotation (Placement(
            transformation(extent={{-58,18},{18,48}})),
          __Dymola_choicesAllMatching=true);
      replaceable
      Components.AdanVenousRed._7af7a4.HeartComponent heartComponent(
        UseFrequencyInput=true,
        UseThoracicPressureInput=true,
        HR=1)
        constrainedby Components.Auxiliary.HeartBase
        annotation (Placement(transformation(extent={{-16,-32},{-36,-12}})));
      replaceable
      Components.AdanVenousRed.PulmonaryComponent pulmonaryComponent(
          UseThoracic_PressureInput=true,                            pulmonary(
          u_pas(start=3871.5508),
          u_pat(start=3871.314),
          u_par(start=3863.3025),
          u_pcp(start=3634.0552),
          u_pvn(start=1266.1965),
          v_pas(start=1.0003076e-06),
          v_pat(start=2.2090626e-05))) constrainedby
        Components.Auxiliary.PulmonaryBase
        annotation (Placement(transformation(extent={{-34,-62},{-14,-42}})));
    Modelica.Blocks.Sources.Trapezoid phi(
        amplitude=0.75,
        rising=2,
        width=20,
        falling=2,
        period=200,
        nperiod=1,
        offset=0.25,
        startTime=20)
        annotation (Placement(transformation(extent={{-70,-6},{-50,14}})));
      inner Components.Settings settings
        annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      Components.ConditionalConnection condHR(disconnectedValue=0.25, disconnected=true)
                             annotation (Placement(transformation(extent={{42,55.2592},
                {54,65.9259}})));
      Components.Baroreflex.HeartRate heartRate
        annotation (Placement(transformation(extent={{16,-28},{4,-16}})));
    Modelica.Blocks.Sources.Trapezoid           thoracic_pressure(
        amplitude=40*133,
        rising=2,
        width=20,
        falling=2,
        period=200,
        nperiod=1,
        offset=0,
        startTime=40)
        annotation (Placement(transformation(extent={{-102,-50},{-82,-30}})));
      Components.ConditionalConnection condTP(disconnectedValue=0, disconnected=true)
                 annotation (Placement(transformation(extent={{-69,-45.3333},{-55,-33.3333}})));
      Components.ConditionalConnection condPhi(disconnectedValue=0.25, disconnected=
           false)             annotation (Placement(transformation(extent={{-26,-0.44444},
                {-14,9.5556}})));
    equation
      connect(Systemic1.port_b, heartComponent.sv) annotation (Line(
          points={{18,28},{38,28},{38,-12},{-16,-12}},
          color={0,0,0},
          thickness=1));
      connect(heartComponent.sa, Systemic1.port_a) annotation (Line(
          points={{-36,-12},{-76,-12},{-76,28},{-58,28}},
          color={0,0,0},
          thickness=1));
      connect(heartComponent.pv, pulmonaryComponent.port_b) annotation (Line(
          points={{-16,-32},{6,-32},{6,-52},{-14,-52}},
          color={0,0,0},
          thickness=1));
      connect(pulmonaryComponent.port_a, heartComponent.pa) annotation (Line(
          points={{-34,-52},{-52,-52},{-52,-32},{-36,-32}},
          color={0,0,0},
          thickness=1));
      connect(condHR.y,heartRate. phi) annotation (Line(points={{54.6,60},{68,
              60},{68,-22},{16,-22}},color={0,0,127}));
      connect(heartRate.HR, heartComponent.frequency_input)
        annotation (Line(points={{3.88,-22},{-16,-22}},  color={0,0,127}));
      connect(thoracic_pressure.y,condTP. u) annotation (Line(points={{-81,-40},
              {-76,-40},{-76,-40},{-70.4,-40}}, color={0,0,127}));
      connect(condTP.y, Systemic1.thoracic_pressure_input) annotation (Line(
            points={{-54.3,-40},{-46,-40},{-46,20},{-28,20}}, color={0,0,127}));
      connect(condTP.y, heartComponent.thoracic_pressure_input) annotation (
          Line(points={{-54.3,-40},{-26,-40},{-26,-32}}, color={0,0,127}));
      connect(condTP.y, pulmonaryComponent.thoracic_pressure) annotation (Line(
            points={{-54.3,-40},{-46,-40},{-46,-62},{-24,-62}}, color={0,0,127}));
      connect(Systemic1.phi_input,condPhi. y) annotation (Line(points={{-14,20},{-4,
              20},{-4,4.00002},{-13.4,4.00002}},
                                      color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=60,
          Interval=0.02,
          __Dymola_Algorithm="Cvode"));
    end CVS_7af;
  annotation(preferredView="info",
  version="2.3.2-beta",
  versionBuild=1,
  versionDate="2015-09-15",
  dateModified = "2015-09-15 12:49:00Z",
  revisionId="",
  uses(Modelica(version="3.2.2")), Documentation(revisions="<html>
<p>Generated from <a href=\"https://models.cellml.org/workspace/4ac\">https://models.cellml.org/workspace/4ac</a> Revision: b580e909bfa88dbf598e9fd1f4b15024e676e9b6 from Date: 2019-04-08 8:13:36 AM, Message: tuning the param for veins</p>
</html>"));
  end AdanVenousRed_Safaei;

  model ADAN_main_AdanVenousRed_Safaei_Experiments_simple_phi_base
   extends ADAN_main.AdanVenousRed_Safaei.Experiments.simple_phi_base;
    annotation(experiment(
    StopTime=400,
    Interval=0.11,
    Tolerance=1e-05,
    __Dymola_Algorithm="Cvode"));
  end ADAN_main_AdanVenousRed_Safaei_Experiments_simple_phi_base;
  annotation (uses(Physiolibrary(version="2.3.2-beta"), Modelica(version=
            "3.2.3")), experiment(
      StopTime=80,
      __Dymola_NumberOfIntervals=1500,
      Tolerance=0.0005,
      __Dymola_Algorithm="Dassl"));
end ADAN_main;
