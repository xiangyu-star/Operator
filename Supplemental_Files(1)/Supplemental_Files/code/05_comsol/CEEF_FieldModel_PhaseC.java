/*
 * CEEF Phase C — Independent Validation Field Model
 * COMSOL 6.4 Java API
 *
 * Phase C: Shows 4 validation scenarios corresponding to the 6 statistical tests
 * using Liu2019 independent chromatin data.
 *
 * 4 scenarios (matching Phase A/B panel layout):
 *   A: Morula-specific coupling (Test1) — only morula has acc-meth coupling
 *      -> sigma_M very tight, hot spot only at morula
 *   B: ZGA-Reset coupling (Test2) — ZGA predicts morula reset
 *      -> F_ZGA + F_entry combined field
 *   C: Acc-correction correlation (Test3) — acc_morula correlates with correction
 *      -> F_entry scaled by acc_morula (gamma_acc = 1.1619 * rho_boost)
 *   D: Cross-validation consistency (Test4/5) — 99.5% consistency
 *      -> Full control field (most robust scenario)
 *
 * Run:
 *   comsolcompile CEEF_FieldModel_PhaseC.java
 *   comsolbatch -inputfile CEEF_FieldModel_PhaseC.class
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;

public class CEEF_FieldModel_PhaseC {

    static final String OUT_DIR = "C:\\comsolfigs_c";
    static final String MPH_OUT = "C:\\comsolfigs_c\\CEEF_FieldModel_PhaseC.mph";

    static final String[] LABELS = {
        "PhaseC_A_test1_morula_specific",
        "PhaseC_B_test2_ZGA_reset",
        "PhaseC_C_test3_acc_correction",
        "PhaseC_D_test456_crossval"
    };
    static final String[] TITLES = {
        "Test1: Morula-Specific Coupling (F_exit only, hot spot at Morula)",
        "Test2: ZGA-Reset Coupling (F_ZGA only, hot spot at 4-cell)",
        "Test3: Acc-Correction (F_entry only, hot spot at 8-cell)",
        "Test4/5/6: Cross-Validation (Full control, all hot spots)"
    };
    // Colormaps — same family as Phase A/B
    static final String[] CMAPS = {"YlOrRd", "viridis", "hot", "plasma"};

    // Phase C validation parameters
    // Test1: morula-specific — tight sigma_M to isolate morula
    static final double SIGMA_M_TIGHT = 0.10;   // tighter than Phase B (0.2132)
    // Test2: ZGA-Reset — ZGA + entry combined
    static final double GAMMA_ZGA_C   = 1.0;
    static final double GAMMA_ACC_C   = 1.1619;
    // Test3: acc-correction — entry boosted by rho
    static final double GAMMA_ACC_BOOST = 1.1619 * (1.0 + 0.18); // scaled by rho=0.18
    // Test4: full control
    static final double GAMMA_ACC_FULL  = 1.1619;
    static final double GAMMA_CLO_FULL  = 1.8940;

    public static void main(String[] args) throws Exception {
        new File(OUT_DIR).mkdirs();
        Model model = buildModel();
        System.out.println("Saving model ...");
        model.save(MPH_OUT);
        System.out.println("Saved: " + MPH_OUT);
        exportFigures(model);
        ModelUtil.disconnect();
        System.out.println("Done. Figures in: " + OUT_DIR);
    }

    static Model buildModel() throws Exception {
        Model model = ModelUtil.create("CEEF_FieldModel_PhaseC");
        model.modelPath("E:\\progress_comsol_analysis");
        model.label("CEEF Phase C — Independent Validation Field");

        // ── Parameters — correction-vector space (same as Phase A/B) ─────────
        model.param().set("z1_oocyte",  "0.9656");
        model.param().set("z2_oocyte",  "0.2348");
        model.param().set("z1_4cell",   "0.5680");
        model.param().set("z2_4cell",   "0.5554");
        model.param().set("z1_8cell",   "0.2150");
        model.param().set("z2_8cell",   "-0.1490");
        model.param().set("z1_morula",  "3.6330");
        model.param().set("z2_morula",  "-0.2790");
        model.param().set("z1_blast",   "0.6060");
        model.param().set("z2_blast",   "1.5399");
        model.param().set("z1_Katt",    "0.59");
        model.param().set("z2_Katt",    "0.04");

        // Physics parameters — same lambda_K as Phase B
        model.param().set("lambda_K",       "0.05");
        model.param().set("sigma_M",        "0.2132");
        model.param().set("sigma_8",        "0.7415");
        model.param().set("sigma_4",        "0.30");
        model.param().set("gamma_acc",      "1.1619");
        model.param().set("gamma_clo",      "1.8940");
        model.param().set("gamma_zga",      "1.0");
        model.param().set("gamma_re",       "1.0");
        model.param().set("gamma_de",       "1.0");

        // Scenario switches
        model.param().set("sw_zga",    "0.0");
        model.param().set("sw_entry",  "1.0");
        model.param().set("sw_exit",   "1.0");
        model.param().set("sigma_M_use", "0.2132");  // can be overridden per scenario

        // ── Component & 2D geometry ───────────────────────────────────────────
        model.modelNode().create("comp1", true);
        model.geom().create("geom1", 2);
        model.geom("geom1").create("r1", "Rectangle");
        model.geom("geom1").feature("r1").set("pos",  new double[]{-0.6, -2.2});
        model.geom("geom1").feature("r1").set("size", new double[]{6.4, 4.4});
        model.geom("geom1").run("fin");

        // ── Mesh ──────────────────────────────────────────────────────────────
        model.mesh().create("mesh1", "geom1");
        model.mesh("mesh1").create("ftri1", "FreeTri");
        model.mesh("mesh1").feature("size").set("hauto", 2);
        model.mesh("mesh1").run();

        // ── Variables ─────────────────────────────────────────────────────────
        model.variable().create("var1");
        model.variable("var1").model("comp1");

        model.variable("var1").set("rM2", "(x-z1_morula)^2+(y-z2_morula)^2");
        model.variable("var1").set("r82", "(x-z1_8cell)^2+(y-z2_8cell)^2");
        model.variable("var1").set("r42", "(x-z1_4cell)^2+(y-z2_4cell)^2");

        // Gaussian envelopes — sigma_M_use can be changed per scenario
        model.variable("var1").set("gM", "exp(-rM2/(2*sigma_M_use^2))");
        model.variable("var1").set("g8", "exp(-r82/(2*sigma_8^2))");
        model.variable("var1").set("g4", "exp(-r42/(2*sigma_4^2))");

        // F_K baseline
        model.variable("var1").set("FK1", "-lambda_K*(x-z1_Katt)");
        model.variable("var1").set("FK2", "-lambda_K*(y-z2_Katt)");

        // Unit vectors
        model.variable("var1").set("d_zga",
            "sqrt((z1_8cell-z1_4cell)^2+(z2_8cell-z2_4cell)^2)");
        model.variable("var1").set("vzga1", "(z1_8cell-z1_4cell)/d_zga");
        model.variable("var1").set("vzga2", "(z2_8cell-z2_4cell)/d_zga");

        model.variable("var1").set("d_entry",
            "sqrt((z1_morula-z1_8cell)^2+(z2_morula-z2_8cell)^2)");
        model.variable("var1").set("ventry1", "(z1_morula-z1_8cell)/d_entry");
        model.variable("var1").set("ventry2", "(z2_morula-z2_8cell)/d_entry");

        model.variable("var1").set("d_exit",
            "sqrt((z1_blast-z1_morula)^2+(z2_blast-z2_morula)^2)");
        model.variable("var1").set("vexit1", "(z1_blast-z1_morula)/d_exit");
        model.variable("var1").set("vexit2", "(z2_blast-z2_morula)/d_exit");

        // Component fields
        model.variable("var1").set("Fzga1",   "sw_zga*gamma_zga*g4*vzga1");
        model.variable("var1").set("Fzga2",   "sw_zga*gamma_zga*g4*vzga2");
        model.variable("var1").set("Fentry1", "sw_entry*gamma_acc*g8*ventry1");
        model.variable("var1").set("Fentry2", "sw_entry*gamma_acc*g8*ventry2");
        model.variable("var1").set("Fexit1",  "sw_exit*(gamma_re+gamma_de)*gM*vexit1");
        model.variable("var1").set("Fexit2",  "sw_exit*(gamma_re+gamma_de)*gM*vexit2");

        // Total field
        model.variable("var1").set("Fu1", "FK1+Fzga1+Fentry1+Fexit1");
        model.variable("var1").set("Fu2", "FK2+Fzga2+Fentry2+Fexit2");
        model.variable("var1").set("Fmag", "sqrt(Fu1^2+Fu2^2)");

        // Active magnitudes for pseudocolor
        model.variable("var1").set("Fmorula_mag", "sqrt(Fexit1^2+Fexit2^2)");
        model.variable("var1").set("Fzga_mag",    "sqrt(Fzga1^2+Fzga2^2)");
        model.variable("var1").set("Fentry_mag",  "sqrt(Fentry1^2+Fentry2^2)");
        model.variable("var1").set("Fcontrol_mag",
            "sqrt((Fentry1+Fexit1)^2+(Fentry2+Fexit2)^2)");

        model.variable("var1").set("Fu1n", "Fu1/(Fmag+1e-9)");
        model.variable("var1").set("Fu2n", "Fu2/(Fmag+1e-9)");

        // ── Physics ───────────────────────────────────────────────────────────
        model.physics().create("c", "CoefficientFormPDE", "geom1");
        model.physics("c").feature("cfeq1").set("c", "1");
        model.physics("c").feature("cfeq1").set("f", "0");
        model.physics("c").feature("cfeq1").set("a", "0");
        model.physics("c").create("dir1", "DirichletBoundary", 1);
        model.physics("c").feature("dir1").selection().all();
        model.physics("c").feature("dir1").set("r", "0");

        // ── Study ─────────────────────────────────────────────────────────────
        model.study().create("std1");
        model.study("std1").label("Stationary");
        model.study("std1").create("stat", "Stationary");
        model.study("std1").run();

        // ── Plot groups & exports ─────────────────────────────────────────────
        buildPlotGroups(model);
        buildExportFeatures(model);

        return model;
    }

    // Color expressions per scenario — each shows the ACTIVE operator's hot spot
    static final String[] COLOR_EXPRS = {
        "Fmorula_mag",  // A: Test1 — exit only → hot spot at MORULA (right)
        "Fzga_mag",     // B: Test2 — ZGA only → hot spot at 4-CELL (upper left)
        "Fentry_mag",   // C: Test3 — entry only → hot spot at 8-CELL (left)
        "Fcontrol_mag"  // D: Test4/5/6 — full control → ALL hot spots
    };

    static void buildPlotGroups(Model model) throws Exception {
        for (int i = 0; i < 4; i++) {
            String pg = "pg" + (i+1);
            model.result().create(pg, "PlotGroup2D");
            model.result(pg).label(TITLES[i]);
            model.result(pg).set("titletype", "manual");
            model.result(pg).set("title", TITLES[i]);
            model.result(pg).set("xlabel", "e1 (entry correction)");
            model.result(pg).set("ylabel", "e2 (exit correction)");

            model.result(pg).create("surf1", "Surface");
            model.result(pg).feature("surf1").set("expr", COLOR_EXPRS[i]);
            model.result(pg).feature("surf1").set("colortable", CMAPS[i]);
            model.result(pg).feature("surf1").set("colorscalemode", "linear");
            model.result(pg).feature("surf1").set("rangecoloractive", "on");
            model.result(pg).feature("surf1").set("rangecolormin", "0");
            model.result(pg).feature("surf1").set("rangecolormax", "0.5");
            model.result(pg).feature("surf1").set("smooth", "internal");

            model.result(pg).create("str1", "Streamline");
            model.result(pg).feature("str1").set("expr", new String[]{"Fu1", "Fu2"});
            model.result(pg).feature("str1").set("posmethod", "start");
            model.result(pg).feature("str1").set("arrowcount", "200");
            model.result(pg).feature("str1").set("color", "white");
            model.result(pg).feature("str1").set("linetype", "tube");
            model.result(pg).feature("str1").set("tuberadiusscaleactive", "on");
            model.result(pg).feature("str1").set("tuberadiusscale", "0.035");
            model.result(pg).feature("str1").set("inheritcolor", "off");

            model.result(pg).create("arw1", "ArrowSurface");
            model.result(pg).feature("arw1").set("expr", new String[]{"Fu1n", "Fu2n"});
            model.result(pg).feature("arw1").set("xnumber", "22");
            model.result(pg).feature("arw1").set("ynumber", "18");
            model.result(pg).feature("arw1").set("arrowlength", "normalized");
            model.result(pg).feature("arw1").set("arrowscale", "0.30");
            model.result(pg).feature("arw1").set("color", "gray");
            model.result(pg).feature("arw1").set("inheritcolor", "off");
        }
    }

    static void buildExportFeatures(Model model) throws Exception {
        for (int i = 0; i < 4; i++) {
            String expTag = "exp" + (i+1);
            String pg     = "pg" + (i+1);
            String pngPath = OUT_DIR + "\\" + LABELS[i] + ".png";
            model.result().export().create(expTag, pg, "Image2D");
            model.result().export(expTag).set("sourceobject", pg);
            model.result().export(expTag).set("pngfilename", pngPath);
            model.result().export(expTag).set("imagetype",   "png");
            model.result().export(expTag).set("size",        "current");
            model.result().export(expTag).set("width",       "2400");
            model.result().export(expTag).set("height",      "1800");
            model.result().export(expTag).set("unit",        "px");
            model.result().export(expTag).set("resolution",  "300");
            model.result().export(expTag).set("antialias",   "on");
            model.result().export(expTag).set("zoomextents", "on");
            model.result().export(expTag).set("fontsize",    "22");
            model.result().export(expTag).set("logo2d",      "off");
            model.result().export(expTag).set("axes2d",      "on");
            model.result().export(expTag).set("legend2d",    "on");
            model.result().export(expTag).set("title2d",     "on");
        }
    }

    // Scenario parameters: [sw_zga, sw_entry, sw_exit, sigma_M_use, gamma_acc, gamma_clo]
    static final double[][] SCENARIO_PARAMS = {
        {0.0, 0.0, 1.0, 0.2132, 1.1619, 1.8940},  // A: exit only  → morula hot spot (Test1)
        {1.0, 0.0, 0.0, 0.2132, 1.1619, 1.8940},  // B: ZGA only   → 4-cell hot spot (Test2)
        {0.0, 1.0, 0.0, 0.2132, 1.1619, 1.8940},  // C: entry only → 8-cell hot spot (Test3)
        {0.0, 1.0, 1.0, 0.2132, 1.1619, 1.8940},  // D: full control → all hot spots (Test4/5/6)
    };

    static void exportFigures(Model model) throws Exception {
        System.out.println("Exporting Phase C figures ...");
        for (int i = 0; i < 4; i++) {
            double sw_zga   = SCENARIO_PARAMS[i][0];
            double sw_entry = SCENARIO_PARAMS[i][1];
            double sw_exit  = SCENARIO_PARAMS[i][2];
            double sig_m    = SCENARIO_PARAMS[i][3];
            double g_acc    = SCENARIO_PARAMS[i][4];
            double g_clo    = SCENARIO_PARAMS[i][5];

            model.param().set("sw_zga",      String.valueOf(sw_zga));
            model.param().set("sw_entry",    String.valueOf(sw_entry));
            model.param().set("sw_exit",     String.valueOf(sw_exit));
            model.param().set("sigma_M_use", String.valueOf(sig_m));
            model.param().set("gamma_acc",   String.valueOf(g_acc));
            model.param().set("gamma_clo",   String.valueOf(g_clo));

            model.study("std1").run();

            String pg     = "pg" + (i+1);
            String expTag = "exp" + (i+1);
            try {
                System.out.println("  Plotting " + LABELS[i] + " ...");
                model.result(pg).run();
                System.out.println("  Exporting ...");
                model.result().export(expTag).run();
                System.out.println("  OK: " + OUT_DIR + "\\" + LABELS[i] + ".png");
            } catch (Exception e) {
                System.out.println("  WARN [" + LABELS[i] + "]: " + e.getMessage());
            }
        }
    }
}
