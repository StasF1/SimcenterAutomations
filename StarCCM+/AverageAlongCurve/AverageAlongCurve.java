/**
 *  API:            Simcenter STAR-CCM+ 15.04.010
 *  Project:        https://github.com/StasF1/StarCcmMacros
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

import star.base.neo.*;
import star.base.report.*;
import star.common.*;
import star.meshing.*;
import star.vis.*;

public class AverageAlongCurve extends StarMacro {
  public void execute() {
    Simulation simulation = getActiveSimulation();

    String regionOfPipe = "Assembly 1.big_coll";

    // FIXME: Post-process only requried fields
    String[] fieldNames = new String[] {
      "AbsoluteTotalPressure",
      "AbsolutePressure",
      "Density",
      "DynamicViscosity",
      "MassFlux",
      "Temperature",
      "TotalTemperature"
    };

    String pathToSaveCsv = MkdirFromSimulationName(simulation.getSessionPath(), ".PipeCuts");

    {
      List<double[]> origins =
        ReadNumericCsv(simulation.getSessionDirFile() + "\\TEST_DATA.csv");
      Multiply(origins, 0.001); // Convert to metres
      // Add(origins, new double[] {0.0, 0.0, 0.0}); // Move coordinate system

      List<double[]> orientations = Difference(origins);

      for (String field : fieldNames) {
        String pipeCutsCsv = ConvertPipeCutsToCsv(
          CreatePipeCuts(simulation, regionOfPipe, field, origins, orientations)
        );
        SaveTextToFile(pathToSaveCsv + "\\" + field + ".csv", pipeCutsCsv);

        simulation.println(field + " CSV field by the tube length:\n" + pipeCutsCsv);
      }
      simulation.println("End");
    }
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

  // TODO: Create a template for any array type
  private static void Print(Simulation simulation, double[] array) {
    boolean isFirst = true;
    simulation.print("[");
    for (double element : array) {
      if (isFirst) {
        simulation.print(element);
        isFirst = false;
      } else {
        simulation.print(",");
        simulation.print(element);
      }
    }
    simulation.print("]");
  }

  // TODO: Create a template for a list with any array type
  private static void Print(Simulation simulation, List<double[]> container) {
    simulation.print("{");
    for (int i = 0; i < container.size(); i++) {
      Print(simulation, container.get(i));
      if (i + 1 != container.size()) {
        simulation.print(",");
      }
    }
    simulation.print("}\n");
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


  /* -------------------------------- List & arrays processing --------------------------------- */

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

  private static double[] Substract(double[] array, double subtrahend) {
    return Add(array, -subtrahend);
  }
  private static double[] Substract(double[] array, double[] subtrahend) {
    for (int i = 0; i < array.length; i++) {
      array[i] -= subtrahend[i];
    }
    return array;
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

  private static List<double[]> Add(List<double[]> container, double[] subtrahend) {
    for (int i = 0; i < container.size(); i++) {
      Add(container.get(i), subtrahend);
    }
    return container;
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
      diff.add(Substract(container.get(i + 1).clone(), container.get(i)));
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

  private static void CreatePlaneSection(Simulation simulation,
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

      simulation.println("Created " + presentationName
                         + " plane section in the region " + regionName);
    }
  }

  private void EditPlaneSection(Simulation simulation,
                                String regionOfPipe, double[] origin, double[] orientation) {
    CreatePlaneSection(simulation, "alongCurveCut", regionOfPipe);

    PlaneSection planeSection =
      ((PlaneSection) simulation.getPartManager().getObject("alongCurveCut"));

    Units units = ((Units) simulation.getUnitsManager().getObject("m"));

    planeSection.getOriginCoordinate().setCoordinate(
      units, units, units,
      new DoubleVector(origin)
    );
    planeSection.getOrientationCoordinate().setCoordinate(
      units, units, units,
      new DoubleVector(orientation)
    );

    planeSection.getInputParts().setQuery(null);
    planeSection.getInputParts().setObjects(
      simulation.getRegionManager().getRegion(regionOfPipe)
    );
  }

  private static void CreateSurfaceAverageReport(Simulation simulation,
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

      simulation.println("Created " + presentationName
                         + " surface average report on the plane section " + planeSectionName);
    }
  }

  private double GetReportValue(Simulation simulation, String reportName, String fieldName) {
    CreateSurfaceAverageReport(simulation, reportName, "alongCurveCut");

    AreaAverageReport surfaceAveragePipeCutReport =
      ((AreaAverageReport) simulation.getReportManager().getReport(reportName));

    PrimitiveFieldFunction field =
      ((PrimitiveFieldFunction) simulation.getFieldFunctionManager().getFunction(fieldName));

    surfaceAveragePipeCutReport.setFieldFunction(field);
    return surfaceAveragePipeCutReport.getValue();
  }

  private List<PipeCut> CreatePipeCuts(Simulation simulation,
                                       String regionOfPipe, String fieldName,
                                       List<double[]> origins, List<double[]> orientations) {
    List<PipeCut> pipeCuts = new ArrayList<PipeCut>();

    for (int i = 0; i < orientations.size(); i++) {
      double[] origin = origins.get(i);
      double[] orientation = orientations.get(i);

      EditPlaneSection(simulation, regionOfPipe, origin, orientation);
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
}

/* (C) 2021 Stanislav Stasheuskii */