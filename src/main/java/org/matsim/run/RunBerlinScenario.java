/* *********************************************************************** *
 * project: org.matsim.*
 *                                                                         *
 * *********************************************************************** *
 *                                                                         *
 * copyright       : (C) 2017 by the members listed in the COPYING,        *
 *                   LICENSE and WARRANTY file.                            *
 * email           : info at matsim dot org                                *
 *                                                                         *
 * *********************************************************************** *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *   See also COPYING, LICENSE and WARRANTY file                           *
 *                                                                         *
 * *********************************************************************** */

package org.matsim.run;

import ch.sbb.matsim.routing.pt.raptor.SwissRailRaptorModule;
import com.google.inject.Singleton;
import org.apache.log4j.Logger;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.matsim.analysis.RunPersonTripAnalysis;
import org.matsim.api.core.v01.Id;
import org.matsim.api.core.v01.Scenario;
import org.matsim.api.core.v01.TransportMode;
import org.matsim.api.core.v01.network.Link;
import org.matsim.api.core.v01.population.*;
import org.matsim.contrib.drt.routing.DrtRoute;
import org.matsim.contrib.drt.routing.DrtRouteFactory;
import org.matsim.core.config.Config;
import org.matsim.core.config.ConfigGroup;
import org.matsim.core.config.ConfigUtils;
import org.matsim.core.config.groups.PlanCalcScoreConfigGroup.ActivityParams;
import org.matsim.core.config.groups.PlansCalcRouteConfigGroup;
import org.matsim.core.config.groups.QSimConfigGroup.TrafficDynamics;
import org.matsim.core.config.groups.VspExperimentalConfigGroup;
import org.matsim.core.controler.AbstractModule;
import org.matsim.core.controler.Controler;
import org.matsim.core.controler.OutputDirectoryHierarchy;
import org.matsim.core.controler.OutputDirectoryLogging;
import org.matsim.core.gbl.Gbl;
import org.matsim.core.gbl.MatsimRandom;
import org.matsim.core.population.PopulationUtils;
import org.matsim.core.population.routes.RouteFactories;
import org.matsim.core.router.AnalysisMainModeIdentifier;
import org.matsim.core.router.TripStructureUtils;
import org.matsim.core.scenario.ScenarioUtils;
import org.matsim.core.scoring.functions.ScoringParametersForPerson;
import org.matsim.core.utils.geometry.transformations.TransformationFactory;
import org.matsim.core.utils.gis.ShapeFileReader;
import org.matsim.prepare.population.AssignIncome;
import org.matsim.run.drt.OpenBerlinIntermodalPtDrtRouterModeIdentifier;
import org.matsim.run.drt.RunDrtOpenBerlinScenario;
import playground.vsp.scoring.IncomeDependentUtilityOfMoneyPersonScoringParameters;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.matsim.core.config.groups.ControlerConfigGroup.RoutingAlgorithmType.FastAStarLandmarks;

/**
* @author ikaddoura
*/

public final class RunBerlinScenario {

	private static final Logger log = Logger.getLogger(RunBerlinScenario.class );

	public static void main(String[] args) {
		
		for (String arg : args) {
			log.info( arg );
		}
		
		if ( args.length==0 ) {
			args = new String[] {"scenarios/berlin-v5.5-1pct/input/berlin-v5.5-1pct.config.xml"}  ;
		}

		Config config = prepareConfig( args ) ;
		Scenario scenario = prepareScenario( config ) ;
		Controler controler = prepareControler( scenario ) ;
		controler.run();
	}

	public static Controler prepareControler( Scenario scenario ) {
		// note that for something like signals, and presumably drt, one needs the controler object
		
		Gbl.assertNotNull(scenario);
		
		final Controler controler = new Controler( scenario );
		
		if (controler.getConfig().transit().isUseTransit()) {
			// use the sbb pt raptor router
			controler.addOverridingModule( new AbstractModule() {
				@Override
				public void install() {
					install( new SwissRailRaptorModule() );
				}
			} );
		} else {
			log.warn("Public transit will be teleported and not simulated in the mobsim! "
					+ "This will have a significant effect on pt-related parameters (travel times, modal split, and so on). "
					+ "Should only be used for testing or car-focused studies with a fixed modal split.  ");
		}
		
		
		
		// use the (congested) car travel time for the teleported ride mode
		controler.addOverridingModule( new AbstractModule() {
			@Override
			public void install() {
				addTravelTimeBinding( TransportMode.ride ).to( networkTravelTime() );
				addTravelDisutilityFactoryBinding( TransportMode.ride ).to( carTravelDisutilityFactoryKey() );
				bind(AnalysisMainModeIdentifier.class).to(OpenBerlinIntermodalPtDrtRouterModeIdentifier.class);
				
				//use income-dependent marginal utility of money for scoring
				bind(ScoringParametersForPerson.class).to(IncomeDependentUtilityOfMoneyPersonScoringParameters.class).in(Singleton.class);
			}
		} );

		return controler;
	}
	
	public static Scenario prepareScenario( Config config ) {
		Gbl.assertNotNull( config );
		
		// note that the path for this is different when run from GUI (path of original config) vs.
		// when run from command line/IDE (java root).  :-(    See comment in method.  kai, jul'18
		// yy Does this comment still apply?  kai, jul'19

		/*
		 * We need to set the DrtRouteFactory before loading the scenario. Otherwise DrtRoutes in input plans are loaded
		 * as GenericRouteImpls and will later cause exceptions in DrtRequestCreator. So we do this here, although this
		 * class is also used for runs without drt.
		 */
		final Scenario scenario = ScenarioUtils.createScenario( config );

		RouteFactories routeFactories = scenario.getPopulation().getFactory().getRouteFactories();
		routeFactories.setRouteFactory(DrtRoute.class, new DrtRouteFactory());
		
		ScenarioUtils.loadScenario(scenario);

		BerlinExperimentalConfigGroup berlinCfg = ConfigUtils.addOrGetModule(config, BerlinExperimentalConfigGroup.class);
		if (berlinCfg.getPopulationDownsampleFactor() != 1.0) {
			downsample(scenario.getPopulation().getPersons(), berlinCfg.getPopulationDownsampleFactor());
		}

		AssignIncome.assignIncomeToPersonSubpopulationAccordingToGermanyAverage(scenario.getPopulation());

		/* ######################## */

		var network = scenario.getNetwork();	// read scenario
		//var shapeFileName = "shp/Ring_4326.shp";			// read shapefile
		var shapeFileName = "shp/Ring_31468.shp";			// read shapefile
		ArrayList<String> affectedLinkIds = new ArrayList<>();
		var features = ShapeFileReader.getAllFeatures(shapeFileName);
		var geometries = features.stream()
				.map(simpleFeature -> (Geometry) simpleFeature.getDefaultGeometry())
				.collect(Collectors.toList());
		var geometry = geometries.get(0);		// geometry containing polygon at Str. d. 17. Juni
		GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();
		// iterate over all links of the network to check if they are at Str. d. 17. Juni
		for (Link link : network.getLinks().values()) {
			// intersect polygon with link center point
			if (isCovered(geometry, link, geometryFactory) && link.getAllowedModes().contains("car")) {
				affectedLinkIds.add(link.getId().toString());
				Set<String> modes = new HashSet<>(link.getAllowedModes());
				modes.remove("car");link.setAllowedModes(modes);

			}
		}

		Population population = scenario.getPopulation();

		// ANSATZ 1:

		// Fehler: findet keine Route von Node x zu Node y für mode car

		// Unser Code: iteriert über alle person -> plan -> trip -> legs -> route links
		// Wenn ein leg (car) über einen gesperrten Link führt, werden alle activities und legs des trips gelöscht,
		// bis auf den ersten Leg. Dessen Attribute werden gelöscht und der mode auf "car" gesetzt.
		// Problem dabei vlt.: Geht der Hinweg durch gesperrten S-Bahnring, kann nicht mehr der Modus car
		// genutzt werden, stattdessen dann pt/bike. Angenommen, der Rückweg geht außen um den S-Bahnring und die Legs
		// werden somit nicht gelöscht. Dann ist das Auto für den Rückweg nicht verfügbar, da für den Hinweg jetzt z.B.
		// pt genutzt wurde statt car => Fehler.

		int i = 0;
		int max = population.getPersons().values().size();
		for (Person person : population.getPersons().values()) {
			// for seeing progress...
			i++;
			if (i % 100 == 0) {
				System.out.println("Person " + i + "/" + max + "(" + Math.round((float) i / (float) max * 100) + "%)" + ": " + person.getId().toString());
			}
			for (Plan plan : person.getPlans()){
				List<PlanElement> tripElementsToBeRemoved = new ArrayList<>();
				List<TripStructureUtils.Trip> trips = TripStructureUtils.getTrips(plan);
				for (TripStructureUtils.Trip trip : trips){
					boolean removeTripLegs = false;
					for (Leg leg : trip.getLegsOnly()){
						if (leg.getMode().equals("car") && leg.getRoute() != null) {
							String[] linksUsed = leg.getRoute().getRouteDescription().split(" ");
							for (String linkUsed : linksUsed) {
								if (affectedLinkIds.contains(linkUsed)) {
									removeTripLegs = true;
									break;
								}
							}
						}
						else if (leg.getRoute() == null) {
							System.out.println("Leg has no route! Continuing...");
						}
						if (removeTripLegs){ break; }
					}
					if (removeTripLegs){
						tripElementsToBeRemoved.addAll(trip.getTripElements());
						Leg newLeg = PopulationUtils.createLeg("walk");
						int index = PopulationUtils.getActLegIndex(plan, trip.getTripElements().get(0));
						plan.getPlanElements().add(index, newLeg);
					}
				}
				for (PlanElement tripElement: tripElementsToBeRemoved){
					// Wir haben das auch mit PopulationUtils.removeActivity() und mit einer angepassten Methode
					// removeActivitySetPrevLegCar versucht (s.u.). Ist äquivalent zu dem was wir in Z. 224-227 und
					// Z. 234 gemacht haben
					plan.getPlanElements().remove(PopulationUtils.getActLegIndex(plan, tripElement));
					/*
					if (tripElement.toString().contains("interaction")){
						removeActivitySetPrevLegCar(plan, PopulationUtils.getActLegIndex(plan, tripElement));
					}
					 */
				}
			}
		}

		// 	ANSATZ 2:
		// HABEN WIR VERWORFEN.

		// Kein Fehler, aber ca. 10% weniger Trips im Output als im Base Scenario (161.000 ggü. 182.000)

		// Unser Code: iteriert über alle person -> plan -> leg -> route links
		// Wenn ein leg (car) über einen gesperrten Link führt, werden alle legs des plans(!) gelöscht.
		/*
		int i = 0;
		int max = population.getPersons().values().size();
		for (Person person : population.getPersons().values()) {
			// for seeing progress...
			i++;
			if (i % 100 == 0) {
				System.out.println("Person " + i + "/" + max + "(" + (float) i / (float) max * 100 + "%)" + ": " + person.getId().toString());
			}
			//
			for (Plan plan : person.getPlans()) {
				boolean removePlanLegs = false;
				List<Leg> legs = PopulationUtils.getLegs(plan);
				for (Leg leg : legs) {
					if (leg.getMode().equals("car")) {
						try {
							String routeDescription = leg.getRoute().getRouteDescription();
							String[] linksUsed = routeDescription.split(" ");
							for (String linkUsed : linksUsed) {
								if (!removePlanLegs && affectedLinkIds.contains(linkUsed)) {
									//PopulationUtils.removeLeg(plan, PopulationUtils.getActLegIndex(plan, leg));
									removePlanLegs = true;
									break;
								}
							}
						}
						catch (NullPointerException e){
							// tritt nicht mehr auf.
							System.out.println("NullPointerException: Route for Leg not found! Continuing...");
						}
					}
				}
				if (removePlanLegs){
					for (Leg leg : legs){
						PopulationUtils.removeLeg(plan, PopulationUtils.getActLegIndex(plan, leg));
					}
				}
			}
		}*/
		return scenario;
	}

	public static Config prepareConfig( String [] args, ConfigGroup... customModules ){
		return prepareConfig( RunDrtOpenBerlinScenario.AdditionalInformation.none, args, customModules ) ;
	}
	public static Config prepareConfig( RunDrtOpenBerlinScenario.AdditionalInformation additionalInformation, String [] args,
					    ConfigGroup... customModules ) {
		OutputDirectoryLogging.catchLogEntries();
		
		String[] typedArgs = Arrays.copyOfRange( args, 1, args.length );
		
		ConfigGroup[] customModulesToAdd;
		if (additionalInformation == RunDrtOpenBerlinScenario.AdditionalInformation.acceptUnknownParamsBerlinConfig) {
			customModulesToAdd = new ConfigGroup[]{new BerlinExperimentalConfigGroup(true)};
		} else {
			customModulesToAdd = new ConfigGroup[]{new BerlinExperimentalConfigGroup(false)};
		}
		ConfigGroup[] customModulesAll = new ConfigGroup[customModules.length + customModulesToAdd.length];
		
		int counter = 0;
		for (ConfigGroup customModule : customModules) {
			customModulesAll[counter] = customModule;
			counter++;
		}
		
		for (ConfigGroup customModule : customModulesToAdd) {
			customModulesAll[counter] = customModule;
			counter++;
		}
		
		final Config config = ConfigUtils.loadConfig( args[ 0 ], customModulesAll );
		
		config.controler().setRoutingAlgorithmType( FastAStarLandmarks );
		config.controler().setOverwriteFileSetting(OutputDirectoryHierarchy.OverwriteFileSetting.deleteDirectoryIfExists);

		config.controler().setLastIteration(1);//////////////////////////////////////////////////////////////////////////////


		config.subtourModeChoice().setProbaForRandomSingleTripMode( 0.5 );
		
		config.plansCalcRoute().setRoutingRandomness( 3. );
		config.plansCalcRoute().removeModeRoutingParams(TransportMode.ride);
		config.plansCalcRoute().removeModeRoutingParams(TransportMode.pt);
		config.plansCalcRoute().removeModeRoutingParams(TransportMode.bike);
		config.plansCalcRoute().removeModeRoutingParams("undefined");
		
		config.qsim().setInsertingWaitingVehiclesBeforeDrivingVehicles( true );
				
		// vsp defaults
		config.vspExperimental().setVspDefaultsCheckingLevel( VspExperimentalConfigGroup.VspDefaultsCheckingLevel.info );
		config.plansCalcRoute().setAccessEgressType(PlansCalcRouteConfigGroup.AccessEgressType.accessEgressModeToLink);
		config.qsim().setUsingTravelTimeCheckInTeleportation( true );
		config.qsim().setTrafficDynamics( TrafficDynamics.kinematicWaves );
				
		// activities:
		for ( long ii = 600 ; ii <= 97200; ii+=600 ) {
			config.planCalcScore().addActivityParams( new ActivityParams( "home_" + ii + ".0" ).setTypicalDuration( ii ) );
			config.planCalcScore().addActivityParams( new ActivityParams( "work_" + ii + ".0" ).setTypicalDuration( ii ).setOpeningTime(6. * 3600. ).setClosingTime(20. * 3600. ) );
			config.planCalcScore().addActivityParams( new ActivityParams( "leisure_" + ii + ".0" ).setTypicalDuration( ii ).setOpeningTime(9. * 3600. ).setClosingTime(27. * 3600. ) );
			config.planCalcScore().addActivityParams( new ActivityParams( "shopping_" + ii + ".0" ).setTypicalDuration( ii ).setOpeningTime(8. * 3600. ).setClosingTime(20. * 3600. ) );
			config.planCalcScore().addActivityParams( new ActivityParams( "other_" + ii + ".0" ).setTypicalDuration( ii ) );
		}
		config.planCalcScore().addActivityParams( new ActivityParams( "freight" ).setTypicalDuration( 12.*3600. ) );

		ConfigUtils.applyCommandline( config, typedArgs ) ;

		return config ;
	}
	
	public static void runAnalysis(Controler controler) {
		Config config = controler.getConfig();
		
		String modesString = "";
		for (String mode: config.planCalcScore().getAllModes()) {
			modesString = modesString + mode + ",";
		}
		// remove last ","
		if (modesString.length() < 2) {
			log.error("no valid mode found");
			modesString = null;
		} else {
			modesString = modesString.substring(0, modesString.length() - 1);
		}
		
		String[] args = new String[] {
				config.controler().getOutputDirectory(),
				config.controler().getRunId(),
				"null", // TODO: reference run, hard to automate
				"null", // TODO: reference run, hard to automate
				config.global().getCoordinateSystem(),
				"https://svn.vsp.tu-berlin.de/repos/public-svn/matsim/scenarios/countries/de/berlin/projects/avoev/shp-files/shp-bezirke/bezirke_berlin.shp",
				TransformationFactory.DHDN_GK4,
				"SCHLUESSEL",
				"home",
				"10", // TODO: scaling factor, should be 10 for 10pct scenario and 100 for 1pct scenario
				"null", // visualizationScriptInputDirectory
				modesString
		};
		
		try {
			RunPersonTripAnalysis.main(args);
		} catch (IOException e) {
			log.error(e.getStackTrace());
			throw new RuntimeException(e.getMessage());
		}
	}
	
	private static void downsample( final Map<Id<Person>, ? extends Person> map, final double sample ) {
		final Random rnd = MatsimRandom.getLocalInstance();
		log.warn( "Population downsampled from " + map.size() + " agents." ) ;
		map.values().removeIf( person -> rnd.nextDouble() > sample ) ;
		log.warn( "Population downsampled to " + map.size() + " agents." ) ;
	}

	private static boolean isCovered(Geometry geometry, Link link, GeometryFactory geometryFactory){
		return geometry.contains(
				// MATSims Coord object has to be transformed to a Coordinate object
				geometryFactory.createPoint(
						new Coordinate(
								link.getCoord().getX(),
								link.getCoord().getY()
						)
				)
		);
	}

	private static void removeActivitySetPrevLegCar(Plan plan, int index) {
		if (index % 2 == 0 && index >= 0 && index <= plan.getPlanElements().size() - 1) {
			if (plan.getPlanElements().size() == 1) {
				log.warn("" + plan + "[index=" + index + " only one act. nothing removed]");
			} else if (index == 0) {
				plan.getPlanElements().remove(index + 1);
				plan.getPlanElements().remove(index);
			} else if (index == plan.getPlanElements().size() - 1) {
				plan.getPlanElements().remove(index);
				plan.getPlanElements().remove(index - 1);
			} else {
				Leg prev_leg = (Leg)plan.getPlanElements().get(index - 1);
				prev_leg.setDepartureTimeUndefined();
				prev_leg.setTravelTimeUndefined();
				prev_leg.setRoute((Route)null);
				prev_leg.setMode("car");
				plan.getPlanElements().remove(index + 1);
				plan.getPlanElements().remove(index);
			}
		} else {
			log.warn("" + plan + "[index=" + index + " is wrong. nothing removed]");
		}
	}
}

