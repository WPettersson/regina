
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  iOS User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2016, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

#import "NewDim4TriangulationController.h"
#import "PacketTreeController.h"
#import "ReginaHelper.h"
#import "dim4/dim4exampletriangulation.h"
#import "dim4/dim4triangulation.h"

@interface NewDim4TriangulationController ()
@property (weak, nonatomic) IBOutlet UISegmentedControl *types;
@property (weak, nonatomic) IBOutlet UIView *container;
@property (weak, nonatomic) NewPacketPageViewController *pages;
@end

@implementation NewDim4TriangulationController

- (void)viewDidLoad
{
    self.pages = static_cast<NewPacketPageViewController*>(self.childViewControllers.lastObject);
    [self.pages fillWithPages:@[@"newDim4TriEmpty", @"newDim4TriExample", @"newDim4TriIsosig"]
                 pageSelector:self.types
                   defaultKey:@"NewDim4TriangulationPage"];
}

- (IBAction)create:(id)sender
{
    regina::NPacket* ans = [self.pages create];
    if (ans) {
        self.spec.parent->insertChildLast(ans);
        [self.spec created:ans];
        [self dismissViewControllerAnimated:YES completion:nil];
    }
}

- (IBAction)cancel:(id)sender
{
    [self dismissViewControllerAnimated:YES completion:nil];
}

@end

#pragma mark - Empty page

@implementation NewDim4TriangulationEmptyPage

- (regina::NPacket*)create
{
    regina::NPacket* ans = new regina::Dim4Triangulation();
    ans->setLabel("4-D triangulation");
    return ans;
}

@end

#pragma mark - Example triangulation

typedef regina::Dim4Triangulation* (*Dim4TriangulationCreator)();

/**
 * Represents a single option in the examples picker.
 */
@interface Dim4ExampleTriangulation : NSObject

@property (strong, nonatomic) NSString* name;
@property (assign, nonatomic) Dim4TriangulationCreator creator;

+ (id)exampleWithName:(NSString*)name creator:(Dim4TriangulationCreator)creator;
- (regina::Dim4Triangulation*)create;

@end

@implementation Dim4ExampleTriangulation

+ (id)exampleWithName:(NSString *)name creator:(Dim4TriangulationCreator)creator
{
    Dim4ExampleTriangulation* e = [[Dim4ExampleTriangulation alloc] init];
    if (e) {
        e.name = name;
        e.creator = creator;
    }
    return e;
}

- (regina::Dim4Triangulation *)create
{
    regina::Dim4Triangulation* ans = (*self.creator)();
    ans->setLabel(self.name.UTF8String);
    return ans;
}

@end

#pragma mark - Example page

@interface NewDim4TriangulationExamplePage () <UIPickerViewDataSource, UIPickerViewDelegate> {
    NSArray* options;
}
@property (weak, nonatomic) IBOutlet UIPickerView *example;
@end

#define KEY_LAST_EXAMPLE @"NewDim4TriangulationExample"

@implementation NewDim4TriangulationExamplePage

- (void)viewDidLoad
{
    options = @[[Dim4ExampleTriangulation exampleWithName:@"Minimal 4-sphere (2 pentachora)" creator:&regina::Dim4ExampleTriangulation::fourSphere],
                [Dim4ExampleTriangulation exampleWithName:@"Simplicial 4-sphere (6 pentachora)" creator:&regina::Dim4ExampleTriangulation::simplicialFourSphere],
                [Dim4ExampleTriangulation exampleWithName:@"ℝP⁴" creator:&regina::Dim4ExampleTriangulation::rp4],
                [Dim4ExampleTriangulation exampleWithName:@"Product S³ × S¹" creator:&regina::Dim4ExampleTriangulation::s3xs1],
                [Dim4ExampleTriangulation exampleWithName:@"Twisted product S³ ×~ S¹" creator:&regina::Dim4ExampleTriangulation::s3xs1Twisted],
                [Dim4ExampleTriangulation exampleWithName:@"Cappell-Shaneson knot complement" creator:&regina::Dim4ExampleTriangulation::cappellShaneson]];

    self.example.dataSource = self;
    self.example.delegate = self;
    
    [self.example selectRow:[[NSUserDefaults standardUserDefaults] integerForKey:KEY_LAST_EXAMPLE] inComponent:0 animated:NO];
}

- (NSInteger)numberOfComponentsInPickerView:(UIPickerView *)pickerView
{
    return 1;
}

- (NSInteger)pickerView:(UIPickerView *)pickerView numberOfRowsInComponent:(NSInteger)component
{
    return options.count;
}

- (NSString *)pickerView:(UIPickerView *)pickerView titleForRow:(NSInteger)row forComponent:(NSInteger)component
{
    return [options[row] name];
}

- (void)pickerView:(UIPickerView *)pickerView didSelectRow:(NSInteger)row inComponent:(NSInteger)component
{
    [[NSUserDefaults standardUserDefaults] setInteger:[self.example selectedRowInComponent:0] forKey:KEY_LAST_EXAMPLE];
}

- (regina::NPacket *)create
{
    return [options[[self.example selectedRowInComponent:0]] create];
}

@end

#pragma mark - Isosig page

@interface NewDim4TriangulationIsosigPage ()
@property (weak, nonatomic) IBOutlet UITextField *isosig;
@end

@implementation NewDim4TriangulationIsosigPage

- (IBAction)editingEnded:(id)sender {
    NewDim4TriangulationController* c = static_cast<NewDim4TriangulationController*>(self.parentViewController.parentViewController);
    [c create:sender];
}

- (regina::NPacket *)create
{
    std::string sig = [self.isosig.text stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]].UTF8String;
    if (sig.empty()) {
        UIAlertView* alert = [[UIAlertView alloc]
                              initWithTitle:@"Empty Isomorphism Signature"
                              message:@"Please type an isomorphism signature into the box provided."
                              delegate:nil
                              cancelButtonTitle:@"Close"
                              otherButtonTitles:nil];
        [alert show];
        return 0;
    }
    
    regina::Dim4Triangulation* t = regina::Dim4Triangulation::fromIsoSig(sig);
    if (! t) {
        UIAlertView* alert = [[UIAlertView alloc]
                              initWithTitle:@"Invalid Isomorphism Signature"
                              message:nil
                              delegate:nil
                              cancelButtonTitle:@"Close"
                              otherButtonTitles:nil];
        [alert show];
        return 0;
    }
    
    t->setLabel(sig);
    return t;
}

@end
