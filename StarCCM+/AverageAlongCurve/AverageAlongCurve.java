/**
 *  API:            Simcenter STAR-CCM+ 15.04.010
 *  Project:        https://github.com/StasF1/SimcenterAutomations
 *  License:        GNU General Public License 3.0 ( see LICENSE )
 *  Author:         Stanislav Stashevskii
 *  
 *  Macro:          AverageAlongCurve.java
 *  Description:    Avarage parameters along a curve
 */
package macro;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.DoubleStream;

import star.base.neo.*;
import star.base.report.*;
import star.common.*;
import star.meshing.*;
import star.vis.*;

public class AverageAlongCurve extends StarMacro {
  public void execute() {
    Simulation simulation = getActiveSimulation();

    String side = "L";
    String zone = "D";
    String pipeName = "In";
    int regime = 60;

    // FIXME: Process a vector fields magnitude
    String[] fieldNames = new String[] {
      "AbsoluteTotalPressure"
    };

    String pipeDataPath = (
      simulation.getSessionDirFile() + "\\..\\..\\pipes.data\\"
      + side + "\\" + zone + "\\" + pipeName
    );

    // String pathToSaveCsv = MkdirFromSimulationName(simulation.getSessionPath(), ".data");
    String pathToSaveCsv = pipeDataPath + "\\Loads\\" + String.valueOf(regime);
    {
      File pathToCreate = new File(pathToSaveCsv);
      pathToCreate.mkdir();
    }

    AverageAlongCurve_(simulation, pipeDataPath, pathToSaveCsv, fieldNames);
  }

  private void AverageAlongCurve_(Simulation simulation,
                                  String pipeDataPath, String pathToSaveCsv,
                                  String[] fieldNames) {
    String regionName = "Assembly 1.big_coll";

    List<double[]> origins = ReadNumericCsv(pipeDataPath + "\\origins.csv");
    Multiply(origins, 0.001); // Convert to metres

    List<double[]> normals = Normalize(Difference(origins));

    simulation.println("\nList of origins:");
    Print(simulation, origins);
    simulation.println("\nList of normals:");
    Print(simulation, normals);
    simulation.println("");

    for (String field : fieldNames) {
      String pipeCutsCsv = ConvertPipeCutsToCsv(
        CreatePipeCuts(simulation, regionName, field, origins, normals)
      );
      SaveTextToFile(pathToSaveCsv + "\\" + field + ".csv", pipeCutsCsv);

      simulation.println(field + " CSV field by the tube length:\n" + pipeCutsCsv);
    }
    simulation.println("End");
  }

  /* ----------------------------------- Arrays manipulation ----------------------------------- */

  private static double[] Negative(double[] array) {
    for (int i = 0; i < array.length; i++) {
      array[i] *= -1.0;
    }
    return array;
  }

  private static double Magnitude(double[] array) {
    return Math.sqrt(DoubleStream.of(Multiply(array.clone(), array)).sum());
  }

  private static double[] Normalize(double[] array) {
    return Multiply(array, 1.0/Magnitude(array));
  }

  private static List<double[]> Normalize(List<double[]> container) {
    for (int i = 0; i < container.size(); i++) {
      Normalize(container.get(i));
    }
    return container;
  }

  private static double[] Add(double[] array, double summand) {
    for (int i = 0; i < array.length; i++) {
      array[i] += summand;
    }
    return array;
  }
  private static double[] Add(double[] array, double[] summand) {
    for (int i = 0; i < array.length; i++) {
      array[i] += summand[i];
    }
    return array;
  }

  private static List<double[]> Add(List<double[]> container, double[] subtrahend) {
    for (int i = 0; i < container.size(); i++) {
      Add(container.get(i), subtrahend);
    }
    return container;
  }

  private static double[] Multiply(double[] array, double multiplier) {
    for (int i = 0; i < array.length; i++) {
      array[i] *= multiplier;
    }
    return array;
  }
  private static double[] Multiply(double[] array, double[] multiplier) {
    for (int i = 0; i < array.length; i++) {
      array[i] *= multiplier[i];
    }
    return array; 
  }

  private static double[] DotProduct(double[] a, double[] b) {
    double[] product = new double[3];
    product[0] = a[1]*b[2] - a[2]*b[1];
    product[1] = a[2]*b[0] - a[0]*b[2];
    product[2] = a[0]*b[1] - a[1]*b[0];
    return product;
  }

  private static double[][] CreateRightHandNormals(double[] k) {
    double[] i = new double[k.length];
    double[] j = new double[k.length];
    double[] kDescartes = new double[] {0.0, 0.0, 1.0};

    if (Arrays.equals(k, kDescartes)) {
      i = new double[] {1.0, 0.0, 0.0};
      j = new double[] {0.0, 1.0, 0.0};
    } else {
      i = DotProduct(k, kDescartes);
      j = DotProduct(k, i);
      Normalize(i);
      Normalize(j);
    }

    return new double[][] {i, j, k};
  }

  private static List<double[]> Multiply(List<double[]> container, double multiplier) {
    for (int i = 0; i < container.size(); i++) {
      Multiply(container.get(i), multiplier);
    }
    return container;
  }

  private static List<double[]> Difference(List<double[]> container) {
    List<double[]> diff = new ArrayList<double[]>();
    for (int i = 0; i + 1 < container.size(); i++) {
      diff.add(Add(container.get(i + 1).clone(), Negative(container.get(i).clone())));
    }
    diff.add(diff.get(diff.size() - 1));
    return diff;
  }


  /* ---------------------------------------- StarCCM+ ----------------------------------------- */

  private class PipeCut {
    public PipeCut(double[] origin, double value) {
      coordinates_ = origin;
      value_ = value;
    }

    public double[] getCoordinates() {
      return coordinates_;
    }
    public double getValue() {
      return value_;
    }

    private double[] coordinates_;
    private double value_;
  };

  private CylindricalCoordinateSystem CreateCylindricalCoordinateSystem(Simulation simulation,
                                                                        String presentationName) {
    LabCoordinateSystem labCoordinateSystem = 
      simulation.getCoordinateSystemManager().getLabCoordinateSystem();

    if (!labCoordinateSystem.getLocalCoordinateSystemManager().has(presentationName)) {
      Units units = simulation.getUnitsManager().getPreferredUnits(
        new IntVector(
          new int[] {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        )
      );

      CylindricalCoordinateSystem cylindricalCoordinateSystem = 
        labCoordinateSystem.getLocalCoordinateSystemManager().createLocalCoordinateSystem(
          CylindricalCoordinateSystem.class, "Cylindrical"
        );

      cylindricalCoordinateSystem.getOrigin().setCoordinate(
        units, units, units,
        new DoubleVector(new double[] {0.0, 0.0, 0.0})
      );
      cylindricalCoordinateSystem.getOrigin().setUnits0(units);
      cylindricalCoordinateSystem.getOrigin().setUnits1(units);
      cylindricalCoordinateSystem.getOrigin().setUnits2(units);

      cylindricalCoordinateSystem.getOrigin().setDefinition("");

      cylindricalCoordinateSystem.getOrigin().setValue(
        new DoubleVector(new double[] {0.0, 0.0, 0.0})
      );
      cylindricalCoordinateSystem.setBasis0(
        new DoubleVector(new double[] {1.0, 0.0, 0.0})
      );
      cylindricalCoordinateSystem.setBasis1(
        new DoubleVector(new double[] {0.0, 1.0, 0.0})
      );

      cylindricalCoordinateSystem.setPresentationName(presentationName);
      return cylindricalCoordinateSystem;
    } else {
      return ((CylindricalCoordinateSystem) labCoordinateSystem.getLocalCoordinateSystemManager()
                                                               .getObject(presentationName));
    }
  }

  private void EditCylindricalCoordinateSystem(Simulation simulation,
                                               double[] origin,
                                               double[] normal) {
    CylindricalCoordinateSystem cylindricalCoordinateSystem = CreateCylindricalCoordinateSystem(
      simulation,
      "pipeCylindrical"
    );

    Units units = ((Units) simulation.getUnitsManager().getObject("m"));
    cylindricalCoordinateSystem.getOrigin().setCoordinate(
      units, units, units,
      new DoubleVector(origin)
    );

    double[][] axes = CreateRightHandNormals(normal);
    cylindricalCoordinateSystem.setBasis0(new DoubleVector(axes[0]));
    cylindricalCoordinateSystem.setBasis1(new DoubleVector(axes[1]));
  }

  private ThresholdPart EditPipeThreshold(Simulation simulation,
                                          String regionName,
                                          String presentationName,
                                          double radius) {
    if (!simulation.getPartManager().has(presentationName)) {
      LabCoordinateSystem labCoordinateSystem = 
        simulation.getCoordinateSystemManager().getLabCoordinateSystem();

      CylindricalCoordinateSystem cylindricalCoordinateSystem = 
        ((CylindricalCoordinateSystem) labCoordinateSystem.getLocalCoordinateSystemManager()
                                                          .getObject("pipeCylindrical"));

      Units units = simulation.getUnitsManager().getPreferredUnits(new IntVector(
          new int[] {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
      ));

      Region region = simulation.getRegionManager().getRegion(regionName);

      PrimitiveFieldFunction primitiveFieldFunction = 
        ((PrimitiveFieldFunction) simulation.getFieldFunctionManager().getFunction("Centroid"));

      VectorComponentFieldFunction vectorComponentFieldFunction = 
        (
          (VectorComponentFieldFunction)
          primitiveFieldFunction.getFunctionInCoordinateSystem(cylindricalCoordinateSystem)
                                .getComponentFunction(0)
        );

      ThresholdPart thresholdPart = simulation.getPartManager().createThresholdPart(
        new NeoObjectVector(new Object[] {region}),
        new DoubleVector(new double[] {0.0, radius}),
        units,
        vectorComponentFieldFunction,
        0
      );
      thresholdPart.setPresentationName(presentationName);
      return thresholdPart;
    } else {
      ThresholdPart thresholdPart = 
        ((ThresholdPart) simulation.getPartManager().getObject(presentationName));

      thresholdPart.getRangeQuantities().setArray(new DoubleVector(new double[] {0.0, radius}));
      return thresholdPart;
    }
  }

  private PlaneSection CreatePlaneSection(Simulation simulation,
                                          String presentationName, String regionName) {
    if (!simulation.getPartManager().has(presentationName)) {
      PlaneSection planeSection = (PlaneSection) simulation.getPartManager().createImplicitPart(
        new NeoObjectVector(new Object[] {}),
        new DoubleVector(new double[] {0.0, 0.0, 1.0}),
        new DoubleVector(new double[] {0.0, 0.0, 0.0}),
        0,
        1,
        new DoubleVector(new double[] {0.0})
      );

      Units units = ((Units) simulation.getUnitsManager().getObject("m"));

      planeSection.getOriginCoordinate().setCoordinate(
        units, units, units,
        new DoubleVector(new double[] {0.0, 0.0, 0.0})
      );

      planeSection.getInputParts().setQuery(null);
      Region region = simulation.getRegionManager().getRegion(regionName);
      planeSection.getInputParts().setObjects(region);
      planeSection.setPresentationName(presentationName);

      return planeSection;
    } else {
      return ((PlaneSection) simulation.getPartManager().getObject(presentationName));
    }
  }

  private PlaneSection EditPlaneSection(Simulation simulation,
                                        String regionName,
                                        double[] origin, double[] normal) {
    PlaneSection planeSection = CreatePlaneSection(simulation, "alongCurveCut", regionName);

    Units units = ((Units) simulation.getUnitsManager().getObject("m"));

    planeSection.getOriginCoordinate().setCoordinate(
      units, units, units,
      new DoubleVector(origin)
    );
    planeSection.getOrientationCoordinate().setCoordinate(
      units, units, units,
      new DoubleVector(normal)
    );

    planeSection.getInputParts().setQuery(null);
    planeSection.getInputParts().setObjects(
      simulation.getRegionManager().getRegion(regionName)
    );
    return planeSection;
  }

  private void EditPlaneSection(Simulation simulation,
                                String regionName,
                                ThresholdPart thresholdPart,
                                double[] origin, double[] normal) {
    PlaneSection planeSection = EditPlaneSection(simulation, regionName, origin, normal);

    planeSection.getInputParts().setObjects(thresholdPart);
  }

  private AreaAverageReport CreateSurfaceAverageReport(Simulation simulation,
                                                       String presentationName,
                                                       String planeSectionName) {
    if (!simulation.getReportManager().has(presentationName)) {
      PlaneSection planeSection = 
        ((PlaneSection) simulation.getPartManager().getObject(planeSectionName));

      AreaAverageReport areaAverageReport = 
        simulation.getReportManager().createReport(AreaAverageReport.class);
      areaAverageReport.getParts().setQuery(null);
      areaAverageReport.getParts().setObjects(planeSection);
      areaAverageReport.setPresentationName(presentationName);
      return areaAverageReport;
    } else {
      return ((AreaAverageReport) simulation.getReportManager().getReport(presentationName));
    }
  }

  private double GetReportValue(Simulation simulation, String reportName, String fieldName) {
    AreaAverageReport surfaceAveragePipeCutReport =
      CreateSurfaceAverageReport(simulation, reportName, "alongCurveCut");

    PrimitiveFieldFunction field =
      ((PrimitiveFieldFunction) simulation.getFieldFunctionManager().getFunction(fieldName));

    surfaceAveragePipeCutReport.setFieldFunction(field);
    return surfaceAveragePipeCutReport.getValue();
  }

  private List<PipeCut> CreatePipeCuts(Simulation simulation,
                                       String regionName, String fieldName,
                                       List<double[]> origins, List<double[]> normals) {
    List<PipeCut> pipeCuts = new ArrayList<PipeCut>();

    for (int i = 0; i < normals.size(); i++) {
      double[] origin = origins.get(i);
      double[] normal = normals.get(i);

      EditCylindricalCoordinateSystem(simulation, origin, normal);
      EditPlaneSection(
        simulation,
        regionName, EditPipeThreshold(simulation, regionName, "pipeThreshold", 0.15),
        origin, normal
      );

      pipeCuts.add(new PipeCut(
          origin,
          GetReportValue(simulation, "surfaceAverageAlongCurveCut", fieldName)
      ));
    }
    return pipeCuts;
  }

  private String ConvertPipeCutsToCsv(List<PipeCut> pipeCuts) {
    String pipeCutsCsv = new String("x,y,z,value\n");

    for (int i = 0; i < pipeCuts.size(); i++) {
      PipeCut pipeCut = pipeCuts.get(i);

      for (double coordinate : pipeCut.getCoordinates()) {
        pipeCutsCsv += String.valueOf(coordinate) + ",";
      }
      pipeCutsCsv += pipeCut.getValue() + "\n";
    }
    return pipeCutsCsv;
  }

  /* --------------------------------------- IO helpers ---------------------------------------- */

  private static List<double[]> ReadNumericCsv(String path) {
    return ReadNumericCsv(path, ",", 1);
  }
  private static List<double[]> ReadNumericCsv(String path, String delimeter) {
    return ReadNumericCsv(path, delimeter, 1);
  }
  private static List<double[]> ReadNumericCsv(String path, int headerLine) {
    return ReadNumericCsv(path, ",", headerLine);
  }
  private static List<double[]> ReadNumericCsv(String path, String delimeter, int headerLine) {
    List<double[]> rows = new ArrayList<double[]>();
    try (BufferedReader br = new BufferedReader(new FileReader(path))) {
      String line;
      int lineCounter = 1;
      while ((line = br.readLine()) != null) {
        if (lineCounter > headerLine) {
          rows.add(
            Arrays.stream(line.split(delimeter)).mapToDouble(Double::parseDouble).toArray()
          );
        }
        lineCounter++;
      }
    } catch (IOException exception) {
      System.out.println("Path do not exist --> [" + path + "]");
    }
    return rows;
  }

  // TODO: Create a template for a list with any array type
  private static void Print(Simulation simulation, List<double[]> container) {
    for (int i = 0; i < container.size(); i++) {
        simulation.println(String.valueOf(i) + ": " + Arrays.toString(container.get(i)));
    }
  }

  private static String MkdirFromSimulationName(String sessionPath) {
    return MkdirFromSimulationName(sessionPath, "");
  }
  private static String MkdirFromSimulationName(String sessionPath, String extension) {
    String sessionPathWoExtension = sessionPath.substring(0, sessionPath.lastIndexOf('.'));

    File pathToCreate = new File(sessionPathWoExtension + extension);
    boolean dirCreated = pathToCreate.mkdir();

    return sessionPathWoExtension + extension;
  }

  private static void SaveTextToFile(String filename, String text) {
    // simulation.getSessionDirFile();
    try (PrintWriter out = new PrintWriter(filename)) {
      out.println(text);
    } catch (IOException exception) {
      System.out.println("Path do not exist --> [" + filename + "]");
    }
  }

}

/* (C) 2021 Stanislav Stasheuskii */